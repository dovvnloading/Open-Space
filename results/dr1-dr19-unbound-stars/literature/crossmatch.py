"""Match Gaia IDs as strings and look for changed IDs or historic names by position."""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ARCSEC = np.pi / (180 * 3600)
FIELDS = ["source_id", "catalogue", "literature_id", "release", "name", "match", "separation_arcsec"]


def vectors(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)])


def number(s):
    return float(s) if s and math.isfinite(float(s)) else 0.


def match(candidates, index):
    by_id = defaultdict(list)
    pos, xyz = [], []
    for ref in index:
        if ref["source_id"]:
            by_id[ref["source_id"]].append(ref)
        if not ref["ra_deg"] or not ref["dec_deg"]:
            continue
        ra, dec = float(ref["ra_deg"]), float(ref["dec_deg"])
        if ref["epoch"]:
            # pmra is mu_alpha*cos(delta). All coordinates are ICRS degrees.
            years = 2016 - float(ref["epoch"])
            ra += years * number(ref["pmra"]) / (3.6e6 * math.cos(math.radians(dec)))
            dec += years * number(ref["pmdec"]) / 3.6e6
        pos.append(ref)
        xyz.append((ra,dec))
    unit_vectors = vectors(*zip(*xyz))
    tree = cKDTree(unit_vectors)
    seen = set()
    for r in candidates:
        sid = r["source_id"]
        for ref in by_id[sid]:
            label = "exact_DR3" if ref["release"] in {"DR3","EDR3"} else "identifier_requires_release_check"
            seen.add((sid,ref["catalogue"],ref["source_id"],ref["release"],ref["name"],label,""))
        if not r.get("ra") or not r.get("dec"):
            continue
        uv = vectors([float(r["ra"])],[float(r["dec"])])[0]
        # A broader, deliberately non-identifying flag for records without epoch
        # covers J2000-to-Gaia motion using the candidate's actual proper motion.
        radius = 5 + 25*math.hypot(number(r.get("pmra")), number(r.get("pmdec")))/1000
        for j in tree.query_ball_point(uv, radius*ARCSEC):
            ref = pos[j]
            sep = 2*math.asin(min(1.,np.linalg.norm(unit_vectors[j]-uv)/2))/ARCSEC
            if ref["epoch"] and sep>3:
                continue
            if ref["source_id"]==sid and ref["release"]=="DR3":
                continue
            seen.add((sid,ref["catalogue"],ref["source_id"],ref["release"],ref["name"],
                      "coordinate_review", f"{sep:.3f}"))
    return sorted(seen)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources",nargs="+",type=Path,help="CSV with source_id, ra, dec, pmra and pmdec (case-insensitive)")
    parser.add_argument("--index",type=Path,default=HERE/"index.csv")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    candidates={}
    for source in args.sources:
        with source.open() as f:
            for row in csv.DictReader(f):
                r={k.lower():v for k,v in row.items()}
                candidates[r["source_id"]]=r
    with args.index.open() as f:
        matched=match(candidates.values(),list(csv.DictReader(f)))
    with args.out.open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(FIELDS)
        w.writerows(matched)
    print(f"{len(candidates)} unique input IDs; {len({r[0] for r in matched})} with an ID or positional lead")


if __name__=="__main__":
    main()
