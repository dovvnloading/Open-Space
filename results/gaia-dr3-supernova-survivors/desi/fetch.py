"""Download Swan et al. v1 and Gaia DR3 columns for each exact DR3 designation."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
CAT = "DESI_DR1_catalogue_1.0.csv.gz"
URL = "https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/Tables/" + CAT
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
COLUMNS = ("source_id parallax parallax_error pmra pmdec ruwe visibility_periods_used "
           "duplicated_source astrometric_params_solved nu_eff_used_in_astrometry "
           "pseudocolour ecl_lat phot_g_mean_mag ipd_frac_multi_peak ipd_gof_harmonic_amplitude").split()
BATCH = 220


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fetch_block(job):
    cache, number, ids = job
    path = cache / f"gaia_{number:03}.csv"
    query = ("SELECT " + ",".join(COLUMNS) + " FROM gaiadr3.gaia_source WHERE source_id IN ("
             + ",".join(ids) + ") ORDER BY source_id")
    params = urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=BATCH, QUERY=query))
    for attempt in range(5):
        try:
            raw = path.read_bytes() if path.exists() else urlopen(TAP + "?" + params, timeout=70).read()
            table = csv.DictReader(io.StringIO(raw.decode("utf-8")))
            if table.fieldnames != COLUMNS:
                raise ValueError("unexpected Gaia header")
            rows = list(table)
            got = [row["source_id"] for row in rows]
            if len(got) != len(set(got)) or set(got) != set(ids):
                raise ValueError(f"batch {number} incomplete: {len(got)} of {len(ids)} IDs")
            path.write_bytes(raw)
            return dict(batch=number, rows=len(rows), id_sha256=digest(",".join(ids).encode()),
                        response_sha256=digest(raw), query_sha256=digest(query.encode()))
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def catalogue(cache):
    path = cache / CAT
    if not path.exists():
        path.write_bytes(urlopen(URL, timeout=120).read())
    raw = path.read_bytes()
    if hashlib.md5(raw).hexdigest() != "860b406a8ee1ded7015e87871c7d6e3e":
        raise ValueError("catalogue differs from publisher's v1 checksum")
    with gzip.open(path, "rt") as f:
        return list(csv.DictReader(f))


def main(cache):
    cache.mkdir(parents=True, exist_ok=True)
    cat = catalogue(cache)
    ids = sorted((r["designation"].split()[-1] for r in cat if r["designation"].startswith("Gaia DR3 ")), key=int)
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate DR3 designations")
    jobs = [(cache, i//BATCH, ids[i:i+BATCH]) for i in range(0, len(ids), BATCH)]
    with ThreadPoolExecutor(6) as pool:
        receipts = list(pool.map(fetch_block, jobs))
    receipt = dict(retrieved_utc=datetime.now(timezone.utc).isoformat(), catalogue_url=URL,
                   catalogue_sha256=digest((cache/CAT).read_bytes()), tap=TAP,
                   release="Gaia DR3", maxrec=BATCH, rows=len(cat), dr3_ids=len(ids),
                   batches=receipts)
    (HERE/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    main(p.parse_args().cache)
