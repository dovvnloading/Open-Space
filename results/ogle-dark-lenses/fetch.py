"""Download OGLE-IV Early Warning System event tables and I-band photometry.

Usage: python fetch.py OUTDIR [YEAR ...]

Each season's `lenses.par` and every event's `phot.dat` are stored under
OUTDIR/YEAR/. Existing files are kept, so the command resumes after interruption.
Source: https://www.astrouw.edu.pl/ogle/ogle4/ews/ (Udalski et al. 2015).
"""
import concurrent.futures as cf
import pathlib
import sys
import time
import urllib.request

BASE = "https://www.astrouw.edu.pl/ogle/ogle4/ews"
FIELD_PREFIX = {"BLG": "blg", "GD": "gd", "DG": "dg"}


def get(url, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "open-space-research/0.2"})
            with urllib.request.urlopen(req, timeout=60) as response:
                return response.read()
        except Exception as error:  # network errors are retried
            last = error
            time.sleep(2 + 3 * attempt)
    raise RuntimeError(f"{url}: {last}")


def events(table):
    for line in table.decode("ascii", "replace").splitlines():
        parts = line.split()
        if parts and parts[0][:4].isdigit() and "-" in parts[0]:
            yield parts[0]


def fetch_event(out, name):
    year, kind, number = name.split("-")
    target = out / year / f"{kind.lower()}-{number}.dat"
    if target.exists() and target.stat().st_size > 0:
        return name, "kept"
    data = get(f"{BASE}/{year}/{kind.lower()}-{number}/phot.dat")
    target.write_bytes(data)
    return name, "new"


def main():
    out = pathlib.Path(sys.argv[1])
    years = sys.argv[2:] or ["2022", "2023", "2024", "2025", "2026"]
    names = []
    for year in years:
        (out / year).mkdir(parents=True, exist_ok=True)
        table = get(f"{BASE}/{year}/lenses.par")
        (out / year / "lenses.par").write_bytes(table)
        names += list(events(table))
    print(len(names), "events", flush=True)
    done = failed = 0
    with cf.ThreadPoolExecutor(6) as pool:
        for future in cf.as_completed([pool.submit(fetch_event, out, n) for n in names]):
            try:
                future.result()
                done += 1
            except Exception as error:
                failed += 1
                print("FAILED", error, flush=True)
            if done % 500 == 0:
                print(done, "done", failed, "failed", flush=True)
    print("finished", done, "ok", failed, "failed")


if __name__ == "__main__":
    main()
