"""Retrieve NEOWISE-R single-exposure photometry for a list of positions.

Usage: python fetch.py POSITIONS.csv OUTDIR [--radius ARCSEC] [--batch N] [--workers N]

POSITIONS.csv needs columns `id,ra,dec` (degrees, J2000). Positions are sent to
the IRSA Gator multi-object service in batches. Each batch is stored as a
gzip-compressed IPAC table whose `cntr_01` column is the 1-based row number
within the batch; `batch_NNNNN.ids` lists the corresponding identifiers.
Completed batches are skipped on rerun.

Catalog: NEOWISE-R Single Exposure (L1b) Source Table, `neowiser_p1bs_psd`
(Mainzer et al. 2014; final 2024 release), IRSA doi:10.26131/IRSA144.
"""
import argparse
import concurrent.futures as cf
import csv
import gzip
import pathlib
import time

import requests

URL = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
COLUMNS = ("mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,w1rchi2,w2rchi2,"
           "qual_frame,qi_fact,saa_sep,moon_masked,cc_flags,sso_flg")


def run_batch(outdir, index, rows, radius):
    target = outdir / f"batch_{index:05d}.tbl.gz"
    if target.exists():
        return index, "kept"
    upload = "|   ra         |   dec        |\n|   double     |   double     |\n"
    upload += "".join(f" {float(r['ra']):13.7f}  {float(r['dec']):13.7f} \n" for r in rows)
    for attempt in range(6):
        try:
            response = requests.post(
                URL,
                data={"catalog": "neowiser_p1bs_psd", "spatial": "Upload", "uradius": str(radius),
                      "uradunits": "arcsec", "outfmt": "1", "selcols": COLUMNS},
                files={"filename": ("positions.tbl", upload)}, timeout=3600)
            text = response.text
            if response.status_code == 200 and "\\RowsRetrieved" in text:
                (outdir / f"batch_{index:05d}.ids").write_text("\n".join(r["id"] for r in rows) + "\n")
                temporary = target.with_suffix(".part")
                with gzip.open(temporary, "wt") as handle:
                    handle.write(text)
                temporary.replace(target)
                return index, "new"
            problem = f"status {response.status_code}: {text[:200]!r}"
        except Exception as error:
            problem = repr(error)
        time.sleep(20 + 40 * attempt)
    return index, "FAILED " + problem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("positions")
    parser.add_argument("outdir")
    parser.add_argument("--radius", type=float, default=5.0)
    parser.add_argument("--batch", type=int, default=300)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(args.positions)))
    batches = [rows[i:i + args.batch] for i in range(0, len(rows), args.batch)]
    print(len(rows), "positions in", len(batches), "batches", flush=True)
    start = time.time()
    with cf.ThreadPoolExecutor(args.workers) as pool:
        jobs = [pool.submit(run_batch, outdir, i, b, args.radius) for i, b in enumerate(batches)]
        for n, job in enumerate(cf.as_completed(jobs), 1):
            index, status = job.result()
            print(f"{n}/{len(batches)} batch {index} {status} {time.time() - start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
