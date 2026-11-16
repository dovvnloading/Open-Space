"""Audit official O4b H1/L1 calibration envelopes and vetted spectral lines.

python detector_context.py CACHE OUTPUT_DIRECTORY
Line and calibration checks are descriptive and do not change the search rank.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path

import numpy as np

from acquire import fetch

DCC = "https://dcc.ligo.org/public/"
CAL_NAME = "LIGO_O4b_cal_uncertainty.tgz"
SOURCES = {
    CAL_NAME: ("0202/T2500288/006/", "3d6e7a5ab65421bb2a2175790a98de5b9fd6e0149845255ce45315cfba808450"),
    "O4abH1lines_v2.csv": ("0201/T2500212/002/", "2789e8ce6627b0ce51028ce37fedbfd1f9286357a8b9efddd54c906509c676e1"),
    "O4abL1lines_v2.csv": ("0201/T2500212/002/", "c78dd729ce934625a298a3921aac0e86f68ee9dc0abe127554762774d4f7da30"),
    "O4_injection_params.html": ("0200/T2500198/003/", "3ff5949acd9406026c1d39199875a7343448423a55e7d863a9080f42cbed5e04"),
}


def acquire_context(cache):
    receipts = {}
    for name, (folder, checksum) in SOURCES.items():
        url = DCC + folder + name
        path = cache / name
        fetch(url, path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
            raise ValueError(f"DCC record changed: {name}")
        receipts[name] = {"url": url, "sha256": checksum}
    return receipts


def spectral_bands(cache, detector, gps):
    """Yield active [left,right] Hz and description, using DCC T2500212-v2 formulas."""
    with (cache / f"O4ab{detector}lines_v2.csv").open() as stream:
        rows = csv.reader(stream)
        next(rows)
        for row in rows:
            spacing, kind, offset, first, last, left, right = map(float, row[:7])
            # Some v2 lines embed time tuples inside the comments, including
            # unquoted commas. The first seven numeric fields are unambiguous.
            description = ",".join(row[7:])
            segments = re.findall(r"\((\d+)\s*,\s*(\d+)\)", description)
            if segments and not any(int(a) <= gps < int(b) for a, b in segments):
                continue
            for harmonic in range(int(first), int(last) + 1):
                freq = harmonic * spacing + offset
                scale = freq / (first * spacing + offset) if kind == 2 else 1.0
                yield freq - left * scale, freq + right * scale, row[7].split('[')[0].strip()


def cw_lines(cache, gps):
    """Conservative Doppler allowance incl. Earth and injected binaries, in Hz.

    Uses f0 and secular fdot at the actual GPS, with 1.1e-4 fractional
    allowance for Earth motion and circular-orbit velocity for pulsars 16/17.
    The applicability of each line also depends on the detector's NO_CW bit.
    """
    html = (cache / "O4_injection_params.html").read_text()
    for row in re.findall(r'<tr>(.*?)</tr>', html, re.S):
        fields = re.findall(r'<td>\s*(.*?)\s*</td>', row, re.S)
        if not fields or not fields[0].isdigit():
            continue
        number, ref, f0, fdot = map(float, fields[:4])
        f = f0 + (gps-ref)*fdot
        width = f * (1.1e-4 + (2*np.pi*2.35/83941.2 if fields[17] == "YES" else 0))
        yield f-width, f+width, f"CW pulsar {int(number)}"


def calibration_records(archive, blocks):
    # The DCC instructs users to take the nearest hourly estimate per IFO.
    # Survey intervals can span >1 estimate. Include every estimate nearest
    # any sample second in each block, and treat values independently by IFO.
    pattern = re.compile(r"(H1|L1)_O4b/calibration_uncertainty_(H1|L1)_(\d+)\.txt")
    times = {"H1": [], "L1": []}
    for member in archive:
        match = pattern.fullmatch(member.name)
        if match and match[1] == match[2]:
            times[match[1]].append(int(match[3]))
    for det in times:
        times[det] = np.array(sorted(times[det]))
    summaries = []
    requested = set()
    for block in blocks:
        seconds = block["gps"] + np.arange(4096)
        for det, epochs in times.items():
            right = np.clip(np.searchsorted(epochs, seconds), 1, len(epochs)-1)
            left_t, right_t = epochs[right-1], epochs[right]
            closest = np.where(seconds-left_t <= right_t-seconds, left_t, right_t)
            ids = np.unique(closest)
            requested.update((det,int(i)) for i in ids)
            summaries.append({"block_gps":block["gps"],"role":block["role"],"detector":det,
                              "epochs_gps":ids.tolist(),
                              "largest_epoch_distance_s":int(np.max(np.abs(seconds-closest)))})
    return summaries, requested


def summarize_calibration(cache, blocks):
    path = cache / CAL_NAME
    with tarfile.open(path, "r|gz") as archive:
        summaries, requested = calibration_records(archive, blocks)
    tables = {}
    with tarfile.open(path, "r|gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            key = member.name.split("_")
            identifier = (key[-2],int(key[-1][:-4]))
            if identifier in requested:
                tables[identifier] = np.loadtxt(io.BytesIO(archive.extractfile(member).read()))
    for row in summaries:
        arrays = [tables[row["detector"], t] for t in row["epochs_gps"]]
        output = {}
        for low, high in [(32,896),(150,150),(350,350)]:
            columns = []
            for table in arrays:
                f = table[:,0]
                grid = np.unique(np.r_[low, f[(f>low)&(f<high)], high])
                columns.append(np.array([np.interp(grid, f, table[:,k]) for k in (1,2,3,4,5,6)]))
            mag, phase, lo_m, lo_p, hi_m, hi_p = np.hstack(columns)
            output[f"{low}-{high}Hz"] = {
                "median_magnitude_range":[float(mag.min()),float(mag.max())],
                "median_phase_range_deg":np.rad2deg([phase.min(),phase.max()]).tolist(),
                "magnitude_envelope_16_84":[float(lo_m.min()),float(hi_m.max())],
                "phase_envelope_16_84_deg":np.rad2deg([lo_p.min(),hi_p.max()]).tolist()}
        row["calibration"] = output
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    sources = acquire_context(args.cache)
    selection = json.loads((args.cache/"search_selection.json").read_text())
    summary = {"sources":sources,"calibration":summarize_calibration(args.cache,selection["blocks"])}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"detector_context.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__ == "__main__":
    main()
