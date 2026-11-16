"""Download every LoTSS DR3 component-level source within 3 arcmin of each LoDeSS position."""

import argparse
import csv
import gzip
import io
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table

ENDPOINT = "https://vo.astron.nl/tap/sync"
QUERY = """SELECT t.rowid,s.Source_Name,s.RA,s.DEC,s.Total_flux,s.E_Total_flux,
s.Peak_flux,s.Isl_rms,s.DC_Maj,s.S_Code
FROM TAP_UPLOAD.target AS t JOIN lotss_dr3.main_sources AS s
ON 1=CONTAINS(POINT('ICRS',s.ra,s.dec),CIRCLE('ICRS',t.raj,t.dej,0.05))"""


def fetch(sources, start, output):
    target = output / f"lotss_{start:05d}.csv.gz"
    if target.exists():
        return
    end = min(start + 500, len(sources))
    table = Table(
        [np.arange(start, end), sources["RA"][start:end], sources["DEC"][start:end]],
        names=["rowid", "raj", "dej"],
    )
    with tempfile.TemporaryDirectory() as tmp:
        upload = Path(tmp) / "positions.xml"
        table.write(upload, format="votable")

        def request():
            return subprocess.check_output(
                [
                    "curl",
                    "-sS",
                    "--retry",
                    "3",
                    "--max-time",
                    "150",
                    "-F",
                    "REQUEST=doQuery",
                    "-F",
                    "LANG=ADQL",
                    "-F",
                    "FORMAT=csv",
                    "-F",
                    "MAXREC=200000",
                    "-F",
                    "UPLOAD=target,param:infile",
                    "-F",
                    f"infile=@{upload}",
                    "-F",
                    f"QUERY={QUERY}",
                    ENDPOINT,
                ]
            )

        for attempt in range(5):
            response = request()
            if response.startswith(b"rowid,source_name,ra,dec,"):
                break
            print(f"Retry {start}: {response[-700:]!r}", flush=True)
            time.sleep(2**attempt)
        if (
            not response.startswith(b"rowid,source_name,ra,dec,")
            or response.count(b"\n") >= 200001
        ):
            raise ValueError(
                f"unexpected or truncated TAP result in chunk {start}: {response[-1000:]!r}"
            )
        rows = list(csv.DictReader(io.StringIO(response.decode())))
        if any(not start <= int(row["rowid"]) < end for row in rows):
            raise ValueError("Unexpected upload row identifier")
        part = target.with_suffix(".tmp")
        part.write_bytes(gzip.compress(response, mtime=0))
        part.replace(target)
        print(f"{start}-{end}: {len(rows)} neighbours", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("lodess_fits")
    ap.add_argument("output_dir", type=Path)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sources = Table.read(args.lodess_fits)
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(
            pool.map(
                lambda start: fetch(sources, start, args.output_dir),
                range(0, len(sources), 500),
            )
        )
