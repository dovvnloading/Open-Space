"""Take a lossless Rubin subset from the MPC Isolated Tracklet File (MPC1992)."""
import argparse
import collections
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.request import urlopen

URL = "https://www.minorplanetcenter.net/iau/ITF/itf.txt.gz"
CUTOFF = dt.date(2026, 11, 17)
WINDOW = 29


def digest(path):
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, help="use a previously downloaded full gzip file")
    parser.add_argument("--output", type=Path, default=Path("."))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    headers = {}
    retrieved = dt.datetime.now(dt.timezone.utc).isoformat()
    tmp = tempfile.TemporaryDirectory()
    source = args.file or Path(tmp.name) / "itf.txt.gz"
    if not args.file:
        with urlopen(URL, timeout=150) as response, source.open("wb") as f:
            headers = {k: response.headers.get(k) for k in ["ETag", "Last-Modified", "Content-Length", "Date"]}
            while block := response.read(1024 * 1024):
                f.write(block)
    rows = collections.defaultdict(list)
    monthly = collections.Counter()
    total = 0
    ps1 = []
    ps1_all = []
    ps1_total = 0
    first_day = CUTOFF - dt.timedelta(days=WINDOW)
    with gzip.open(source, "rb") as f:
        for line in f:
            total += 1
            if line[77:80] == b"F51" and line[15:22] in {b"2026 10", b"2026 11"}:
                y, m, d = line[15:32].split()
                day = dt.date(int(y), int(m), int(float(d)))
                if first_day <= day < CUTOFF:
                    ps1_total += 1
                    ps1_all.append(line)
                    try:
                        h, mi, sec = map(float, line[32:44].split())
                        de, am, asec = map(float, line[45:56].split())
                        ra = 15*(h + mi/60 + sec/3600)
                        dec = (1 if line[44:45] == b"+" else -1)*(de + am/60 + asec/3600)
                        if 40 <= ra < 50 and 10 <= dec < 20:
                            ps1.append(line)
                    except ValueError:
                        # Non-position second lines of MPC1992 records.
                        pass
            if line[77:80] != b"X05" or not b"2025" <= line[15:19] <= b"2026":
                continue
            fields = line[15:32].split()
            day = dt.date(int(fields[0]), int(fields[1]), int(float(fields[2])))
            if day < CUTOFF:
                rows[day].append(line)
                monthly[day.strftime("%Y-%m")] += 1
    first = CUTOFF - dt.timedelta(days=WINDOW)
    end = CUTOFF
    primary = sum(len(v) for d, v in rows.items() if first <= d < end)
    # An ITF release can lag the observation dates. Preserve the original count.
    if primary == 0 and rows:
        end = max(rows) + dt.timedelta(days=1)
        first = end - dt.timedelta(days=WINDOW)
    selected = args.output / "x05.mpc"
    with selected.open("wb") as f:
        for day in sorted(rows):
            if first <= day < end:
                for line in rows[day]:
                    f.write(line)
    ps1_file = args.output / "f51.mpc"
    ps1_file.write_bytes(b"".join(ps1))
    all_file = args.output / "f51-all.mpc"
    all_file.write_bytes(b"".join(ps1_all))
    receipt = dict(ps1_all_sha256=digest(all_file), ps1_date_window_all_sky_rows=ps1_total, ps1_patch_rows=len(ps1),
                   ps1_ra_deg=[40,50], ps1_dec_deg=[10,20], ps1_sha256=digest(ps1_file),url=URL, retrieved_utc=retrieved, used_local_copy=bool(args.file),
                   response_headers=headers, source_sha256=digest(source), compressed_bytes=source.stat().st_size,
                   source_lines=total, cutoff_exclusive=str(CUTOFF), primary_29day_rows=primary,
                   x05_calendar_months=dict(sorted(monthly.items())), start=str(first), end_exclusive=str(end),
                   selected_sha256=digest(selected), selected_bytes=selected.stat().st_size,
                   selected_observation_rows_by_utc={str(d):len(rows[d]) for d in sorted(rows) if first <= d < end})
    (args.output / "input.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))

if __name__ == "__main__":
    main()
