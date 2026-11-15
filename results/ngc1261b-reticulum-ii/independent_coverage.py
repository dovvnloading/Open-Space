"""Check all 70 branch stars in Gaia DR3, S5 DR1, and the later Awad catalog (stdlib)."""
import argparse
import csv
import hashlib
import io
import json
import math
import struct
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from reproduce import SOURCES, separation

S5_URL = "https://zenodo.org/api/records/4695135/files/s5_pdr1_light.fits/content"
S5_SHA = "87680e5327c1dc066438489b6b1547952ba2205a9adbccea2f63a2f55d360a61"
GAVO = "https://dc.zah.uni-heidelberg.de/tap/sync"


def fetch(cache, name, url, digest=None):
    file = cache / name
    if file.exists():
        data = file.read_bytes()
    else:
        with urlopen(url, timeout=60) as response:
            data = response.read()
    if digest and hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Unexpected digest for {name}")
    if not file.exists():
        file.write_bytes(data)
    return data


def fits_header(data, offset):
    values = {}
    for i in range(offset, len(data), 80):
        card = data[i:i+80].decode("ascii")
        key = card[:8].strip()
        if key == "END":
            return values, (i+80+2879)//2880*2880
        if card[8:10] == "= ":
            values[key] = card[10:].split("/")[0].strip().strip("'").strip()
    raise ValueError("Incomplete FITS header")


def run(cache):
    cache.mkdir(parents=True, exist_ok=True)
    ibata = fetch(cache, "ibata.txt", *SOURCES["ibata.txt"]).decode()
    sample = [s.split() for s in ibata.splitlines() if s[:19].strip().isdigit() and s[95:97].strip()=="14"]
    ids = [s[0] for s in sample]
    query = "SELECT source_id,radial_velocity,radial_velocity_error FROM gaia.dr3lite WHERE source_id IN (" + ",".join(ids) + ")"
    gaia_url = GAVO + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    gaia = list(csv.DictReader(io.StringIO(fetch(cache,"gavo-dr3.csv",gaia_url).decode())))
    if set(row["source_id"] for row in gaia) != set(ids) or len(gaia)!=len(sample):
        raise ValueError("Incomplete Gaia DR3 coverage")
    data = fetch(cache,"s5_pdr1_light.fits",S5_URL,S5_SHA)
    _, start = fits_header(data,0)
    header, start = fits_header(data,start)
    layout = struct.Struct(">6d19sdqcdq")
    if (header["TTYPE5"], header["TTYPE6"], int(header["NAXIS1"])) != ("ra", "dec", layout.size):
        raise ValueError("Unexpected S5 DR1 binary table")
    nrow = int(header["NAXIS2"])
    coords = [layout.unpack_from(data,start+i*layout.size)[4:6] for i in range(nrow)]
    mins = []
    for s in sample:
        ra, dec = float(s[1]),float(s[2])
        # All candidates within 2 degrees in declination are guaranteed to include
        # every possible six-arcsecond match. Use an exact spherical distance within it.
        mins.append(min((separation(ra,dec,r,d) for r,d in coords if abs(dec-d)<2),default=math.inf))
    awad_url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/693/A69&-out.max=10000"
    tab = fetch(cache,"awad.tsv",awad_url).decode()
    block = next(part for part in tab.split("#Table") if part.startswith("\tJ_A_A_693_A69_ngc1261:"))
    tabrows = [line for line in block.splitlines() if line.split("\t")[0].isdigit()]
    canonical = "\n".join(line.strip() for line in tabrows).encode()
    if hashlib.sha256(canonical).hexdigest() != "61ddfe019c699721500709312e7bf31da2c790ff00b316ac45e4f3e1aacfafb6":
        raise ValueError("VizieR data rows or ordering changed")
    columns = next(line.split("\t") for line in block.splitlines() if line.startswith("GaiaDR3\t"))
    matches = []
    for line in tabrows:
        star = dict(zip(columns, line.split("\t")))
        if star["GaiaDR3"] in ids:
            matches.append({k:star[k].strip() for k in ("GaiaDR3","RV","e_RV","[Fe/H]","e_[Fe/H]","Pmemb")})
    return dict(awad_b1_stars=len(tabrows),awad_b1_exact_matches=matches,atlas_branch_stars=len(sample),gaia_dr3_rows=len(gaia),
                gaia_dr3_nonblank_rvs=sum(bool(row["radial_velocity"]) for row in gaia),
                s5_dr1_rows=nrow,s5_dr1_matched_within_6_arcsec=sum(d<=6 for d in mins))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,required=True)
    print(json.dumps(run(p.parse_args().cache),indent=2))
