"""Check a MeerKAT-predicted sign in separate Parkes single-pulse PSRFITS records."""
import csv
import hashlib
import sqlite3
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import correlate, correlation_lags


def decode(blob):
    # These particular public search-mode snippets contain one primary header
    # and one SUBINT header. Check the full binary layout before reading it.
    ends = [i + 80 for i in range(0, min(16000, len(blob)), 80)
            if blob[i:i+8] == b"END     "][:2]
    assert len(ends) == 2
    header = {blob[i:i+8].strip().decode(): blob[i+10:i+40].split(b"/")[0].strip()
              for i in range(0, ends[1], 80) if blob[i+8:i+10] == b"= "}
    start = ((ends[1] + 2879) // 2880) * 2880
    channels = int(header["NCHAN"])
    samples = int(float(header["NSBLK"]))
    assert (int(header["NPOL"]), int(header["NBITS"]), channels) == (1, 1, 96)
    assert int(header["NAXIS2"]) == 1
    keys = ("TSUBINT", "OFFS_SUB", "DAT_FREQ", "DAT_WTS", "DAT_SCL", "DAT_OFFS", "DATA")
    forms = ("1D", "1D", "96D", "96E", "96E", "96E", f"{samples*12}B")
    for col, (key, form) in enumerate(zip(keys, forms), 1):
        assert header[f"TTYPE{col}"].strip(b"' ").decode() == key
        assert header[f"TFORM{col}"].strip(b"' ").decode() == form
    freqs = np.frombuffer(blob, offset=start+16, dtype=">f8", count=channels)
    offset = start + 16 + channels * 20
    packed = np.frombuffer(blob, dtype="u1", offset=offset, count=samples*12)
    data = np.unpackbits(packed).reshape(samples, channels).T
    order = np.argsort(freqs)
    return data[order], freqs[order], float(header["TBIN"]) * 1000


def profile(blob, dm):
    data, freq, dt_ms = decode(blob)
    assert dt_ms == 0.25
    samples = data.shape[1]
    shifts = np.rint(4.148808e6 * dm * (freq**-2 - freq[-1]**-2) / dt_ms).astype(int)
    # Truncate at the last sample present in every channel after dedispersion.
    n_full = samples - shifts.max()
    if n_full < 65:
        return None
    aligned = np.take_along_axis(data, np.arange(n_full)[None, :] + shifts[:, None], axis=1)
    aligned = aligned.astype(float) - 0.5
    # Smoothing uses the full rolled array as in the exploratory reduction, then
    # accepts only centres at least 32 samples from either coverage edge.
    full = np.take_along_axis(data, (np.arange(samples)[None, :] + shifts[:, None]) % samples, axis=1)
    peak = np.argmax(gaussian_filter1d(full[32:64].mean(0), 3, mode="wrap")[:n_full])
    if peak < 32 or peak >= n_full - 32:
        return None
    return aligned[:, peak-32:peak+33].reshape(3, 32, 65).mean(axis=1)


def lag(a, b):
    vals = correlate(a - a.mean(), b - b.mean())
    lags = correlation_lags(len(a), len(b))
    sel = np.flatnonzero(abs(lags) < 14)
    i = sel[np.argmax(vals[sel])]
    frac = (vals[i-1] - vals[i+1]) / (2 * (vals[i-1] - 2*vals[i] + vals[i+1]))
    return (lags[i] + frac) * .25


def contrast(profiles):
    # Use only the middle third of the receiver to measure brightness. Measure the
    # lag exclusively between the disjoint lower and upper thirds of the band.
    e = profiles[:, 1, 20:45].sum(axis=1)
    cuts = np.quantile(e, [.33, .67])
    faint = profiles[e <= cuts[0]].mean(axis=0)
    bright = profiles[e >= cuts[1]].mean(axis=0)
    return lag(bright[0], bright[2]) - lag(faint[0], faint[2])


def main():
    database, output = map(Path, sys.argv[1:3])
    expected = "12efb081419fb22473b8c0b947cb44c2fbecf2a1b023a9f00b3fded1564fc81a"
    with database.open("rb") as f:
        assert hashlib.file_digest(f, "sha256").hexdigest() == expected
    db = sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True)
    rng = np.random.default_rng(1214)
    rows = []
    for name in ("J1825-0935", "J1840-0815", "J1840-0809"):
        pid, dm = db.execute("SELECT pulsarID,dm FROM pulsar WHERE jname=?", (name,)).fetchone()
        groups = defaultdict(list)
        counts = defaultdict(int)
        for file_id, mult, blob in db.execute(
                "SELECT pfLinkID,pulseNumber,data FROM fileSegment WHERE pulsarID=? ORDER BY segID", (pid,)):
            counts[file_id] += 1
            if mult == 1:
                p = profile(blob, dm)
                if p is not None:
                    groups[file_id].append(p)
        for file_id, total in counts.items():
            filename, mjd = db.execute("SELECT filename,timeStartMJD FROM file WHERE pfLinkID=?",
                                       (file_id,)).fetchone()
            p = np.array(groups[file_id])
            record = dict(pulsar=name, file=filename, mjd=mjd, snippets=total, usable=len(p),
                          bright_minus_faint_lag_ms="", bootstrap_low="", bootstrap_high="")
            if len(p) >= 10:
                draws = [contrast(p[rng.integers(0, len(p), size=len(p))]) for _ in range(1000)]
                lo, hi = np.quantile(draws, [.025, .975])
                record.update(bright_minus_faint_lag_ms=contrast(p), bootstrap_low=lo, bootstrap_high=hi)
            rows.append(record)
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
