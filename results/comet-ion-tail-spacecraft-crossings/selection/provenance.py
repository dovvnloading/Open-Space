"""Record the byte hashes and query URLs of the survey artifacts used in the screen."""

import csv
import hashlib
from datetime import datetime, timedelta

from net import CACHE, ROOT, sburl, craft, hurl


def main():
    rows = []

    def record(name, url):
        p = CACHE / name
        if not p.is_file():
            return
        data = p.read_bytes()
        rows.append(
            dict(
                file=name,
                sha256=hashlib.sha256(data).hexdigest(),
                bytes=len(data),
                url=url,
            )
        )

    record("sbdb.json", sburl)
    for k, (cmd, a, b) in craft.items():
        record(k + ".txt", hurl(cmd, a, b, "1 d"))
    candidates = list(csv.DictReader((CACHE / "screen.csv").open()))
    requests = list(csv.DictReader((CACHE / "refinement-status.csv").open()))
    cs = {(c["craft"], c["comet"], c["date"]): c for c in candidates}
    for item in requests:
        z = cs[item["craft"], item["comet"], item["guess"]]
        t = datetime.fromisoformat(z["date"])
        start = (
            (t - timedelta(days=42 + (33 if z["craft"] == "ulysses" else 16)))
            .date()
            .isoformat()
        )
        stop = (t + timedelta(days=42)).date().isoformat()
        n = f"eph_{z['spkid']}_{start}_{stop}.txt"
        if n not in {r["file"] for r in rows}:
            record(n, hurl(f"DES={z['pdes']};CAP;", start, stop, "6 h"))
    with (ROOT / "survey-inputs.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
