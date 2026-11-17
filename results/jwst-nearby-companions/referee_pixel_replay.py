"""External, clean-runner replay of selected JWST/Spitzer FITS reductions.

NOT an offline unit test. Downloads SHA-pinned files and re-executes published
science code; compares new measurements with frozen outputs. Run in CI, not
from the repository's generic network-disabled verification sandbox.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent


def readcsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def process(folder, *args, timeout=2400):
    cmd = [sys.executable, *map(str, args)]
    print("RUN", folder, cmd, flush=True)
    subprocess.run(cmd, cwd=folder, check=True, timeout=timeout)


def same_rows(original, fresh, keyfields, fields, relative=.05):
    a, b = readcsv(original), readcsv(fresh)
    assert len(a) == len(b), (original, len(a), len(b))
    key = lambda r: tuple(r[k] for k in keyfields)
    orig = {key(r): r for r in a}
    new = {key(r): r for r in b}
    assert orig.keys() == new.keys(), (original, orig.keys() ^ new.keys())
    for k in orig:
        for f in fields:
            expected, got = float(orig[k][f]), float(new[k][f])
            assert math.isclose(expected, got, rel_tol=relative, abs_tol=1e-9), (
                original, k, f, expected, got)
    print("MATCHED", original, "rows", len(a), "fields", fields, flush=True)


def spitzer(tmp):
    folder = ROOT / "g29" / "spitzer"
    frozen = {r["filename"]: r for r in json.loads(
        (folder / "measurements.json").read_text())}
    manifest = json.loads((folder / "spitzer-query.json").read_text())
    cache = folder / "cache"
    cache.mkdir(exist_ok=True)
    for entry in manifest["files"]:
        file = cache / entry["access_url"].split("/")[-1]
        if not file.exists():
            with urlopen(entry["access_url"], timeout=180) as response, file.open("wb") as out:
                while block := response.read(1 << 20):
                    out.write(block)
        with file.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        assert digest == entry["sha256"], (file, digest, entry["sha256"])
        assert file.stat().st_size == entry["bytes"]
        print("FITS SHA256 verified", file.name, flush=True)
    process(folder, "measure.py")
    now = {r["filename"]: r for r in json.loads((folder / "measurements.json").read_text())}
    assert now.keys() == frozen.keys()
    for name in now:
        a, b = frozen[name], now[name]
        assert a["gaia_ties"] == b["gaia_ties"]
        for field, tol in (("fixed_centroid_distance_arcsec", .025),
                           ("source_offset_fixed_vs_comoving_arcsec", .05)):
            if a[field] is None:
                assert b[field] is None
            else:
                assert abs(a[field] - b[field]) <= tol, (name, field, a[field], b[field])
    print("PASS original Spitzer pixel motion", sorted(now), flush=True)


def g29(tmp):
    folder = ROOT / "g29"
    cache, out = tmp / "g29-cache", tmp / "g29-out"
    process(folder, "fetch.py", cache)
    process(folder, "reduce.py", cache, out)
    same_rows(folder / "out" / "contrast_recovery.csv", out / "contrast_recovery.csv",
              ["filter", "radius_arcsec", "contrast"],
              ["coverage_fraction", "threshold_contrast"], relative=.05)
    orig, new = readcsv(folder / "out" / "peaks.csv"), readcsv(out / "peaks.csv")
    assert len(orig) == len(new), (len(orig), len(new))
    print("PASS original G29-38 MIRI image search with six hashed FITS inputs", flush=True)


def nircam(tmp):
    folder = ROOT / "nircam"
    cache1, cache2, out = tmp / "nrc-8826", tmp / "nrc-5835", tmp / "nrc-out"
    process(folder, "nircam_fetch.py", "--inputs", cache1,
            "--manifest", "go8826_inputs.csv")
    process(folder, "nircam_fetch.py", "--inputs", cache2,
            "--manifest", "mcc124_inputs.csv")
    out.mkdir(exist_ok=True, parents=True)
    for name, host in (("l788_34", ".32483977846125917"),
                       ("ross619", ".26697414579879847")):
        process(folder, "nircam_reduce.py", "--field", name, "--inputs",
                cache1, "--out", out, "--host-jy", host)
        same_rows(folder / (name + "_contrast.csv"), out / (name + "_contrast.csv"),
                  ["radius_arcsec"], ["f95_uJy"], relative=.07)
        old, new = readcsv(folder / (name + "_peaks.csv")), readcsv(out / (name + "_peaks.csv"))
        assert len(old) == len(new), (name, len(old), len(new))
        print("PASS", name, "detector-level RDI and injection curves", flush=True)
    # Companion-vs-RDI-reference diagnostics on both fields.
    archived = readcsv(folder / "nircam_trigger_checks.csv")
    process(folder, "nircam_check_sources.py", "--inputs", cache1, "--archival", cache2)
    now = readcsv(folder / "nircam_trigger_checks.csv")
    assert len(archived) == len(now)
    old_ross = [r for r in archived if r["field"] == "ross619"]
    new_ross = [r for r in now if r["field"] == "ross619"]
    assert len(old_ross) == len(new_ross) == 2
    assert all(float(r["archival_rdi_snr"]) < 5 for r in new_ross)
    print("PASS Ross 619 alternate-library rejection", flush=True)


def gj205(tmp):
    folder = ROOT / "gj205"
    cache, out = tmp / "gj205-cache", tmp / "gj205-out"
    for manifest in ("miri_inputs.csv", "miri_f1550_inputs.csv"):
        process(folder, "../nircam/nircam_fetch.py", "--inputs",
                cache, "--manifest", manifest)
    out.mkdir(exist_ok=True, parents=True)
    for band, jy in (("F1065C", "1.263"), ("F1550C", ".603")):
        process(folder, "miri_reduce.py", "--band", band, "--host-jy", jy,
                "--inputs", cache, "--out", out)
        prefix = band.lower()
        same_rows(folder / (prefix + "_miri_contrast.csv"), out / (prefix + "_miri_contrast.csv"),
                  ["radius_arcsec"], ["f95_uJy"], relative=.07)
        old = readcsv(folder / (prefix + "_miri_peaks.csv"))
        new = readcsv(out / (prefix + "_miri_peaks.csv"))
        assert len(old) == len(new), (band, len(old), len(new))
        print("PASS", band, "coronagraphic mask and injection", flush=True)
    process(folder, "miri_source_checks.py", "--inputs", cache, "--out", out)
    print("PASS GJ 205 two-band pixel reduction", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("case", choices=["spitzer", "g29", "nircam", "gj205"])
    p.add_argument("scratch", type=Path)
    a = p.parse_args()
    a.scratch.mkdir(parents=True, exist_ok=True)
    globals()[a.case](a.scratch)
