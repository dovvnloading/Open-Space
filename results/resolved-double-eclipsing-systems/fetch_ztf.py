"""Fetch the ZTF public DR photometry for specified Gaia positions."""
import argparse
import csv
import concurrent.futures
import pathlib
import urllib.parse
import urllib.request


def fetch(row, component, directory):
    name = row[f"source_id{component}"]
    path = directory / (name + ".csv")
    if path.exists():
        return
    ra, dec = row[f"ra{component}"], row[f"dec{component}"]
    query = urllib.parse.urlencode(dict(POS=f"CIRCLE {ra} {dec} 0.00042", FORMAT="csv", BAD_CATFLAGS_MASK="32768", NOBS_MIN=20))
    base = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
    with urllib.request.urlopen(base + query, timeout=120) as response:
        result = response.read()
    if result.lstrip().startswith(b"<"):
        raise ValueError(f"Invalid response for {name}: {result[:140]!r}")
    path.write_bytes(result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pairs")
    p.add_argument("directory",type=pathlib.Path)
    a = p.parse_args()
    a.directory.mkdir(exist_ok=True, parents=True)
    pairs = list(csv.DictReader(open(a.pairs)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        jobs = [pool.submit(fetch,r,i,a.directory) for r in pairs for i in (1,2)]
        for job in jobs:
            job.result()
