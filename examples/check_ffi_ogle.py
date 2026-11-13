"""Compare recurrent FFI sky peaks with position and period in OGLE-IV LMC catalogs."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from open_space.tess_archive import atomic_bytes, download
from open_space.tess_eb_catalog import separation


BASE = "https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/lmc/"
DIRECTORY = {
    "cep": ["cep1O", "cep1O2O", "cep1O2O3O", "cep1O3O", "cep2O", "cep2O3O",
            "cepF", "cepF1O", "cepF1O2O"],
    "rrlyr": ["RRab", "RRc", "RRd", "aRRd"],
    "ecl": ["ecl", "ell"], "acep": ["acepF", "acep1O"],
    "t2cep": ["t2cep"], "dsct": ["dsct"],
}
RADIUS = 30
PERIOD_TOLERANCE = .02


def degrees(value):
    first, minutes, seconds = (float(v) for v in value.split(":"))
    return (-1 if value.startswith("-") else 1) * (abs(first) + minutes/60 + seconds/3600)


def tables(cache):
    records, receipts = {}, []
    for group, classes in DIRECTORY.items():
        for file in ["ident", *classes]:
            rel = f"{group}/{file}.dat"
            path = cache / rel
            raw = path.read_bytes() if path.exists() else download(BASE+rel, 32_000_000)
            if not path.exists():
                atomic_bytes(path, raw)
            n = 0
            for line in raw.decode().splitlines():
                cols = line.split()
                if not cols or not cols[0].startswith("OGLE-LMC-"):
                    continue
                n += 1
                if file == "ident":
                    k = next(j for j,t in enumerate(cols) if t.count(":")==2)
                    records[cols[0]] = {"name": cols[0], "catalog_group": group,
                                        "ra_deg": 15*degrees(cols[k]), "dec_deg": degrees(cols[k+1]),
                                        "pulsation_or_binary_type": " ".join(cols[1:k])}
                else:
                    if cols[0] not in records:
                        continue
                    record = records[cols[0]]
                    # First period appears at column 4. Multi-mode tables use a repeated
                    # block of 8 columns per reported period.
                    periods = [float(cols[3])] if cols[3] != "-" else []
                    if group in {"cep", "rrlyr"} and file in {
                            "cep1O2O", "cep1O2O3O", "cep1O3O", "cep2O3O",
                            "cepF1O", "cepF1O2O", "RRd", "aRRd"}:
                        periods.extend(float(cols[j]) for j in (11,19) if j < len(cols) and cols[j] != "-")
                    record["period_days"] = periods
            receipts.append({"url": BASE+rel, "sha256": sha256(raw).hexdigest(), "rows": n})
    return list(records.values()), receipts


def match(evidence, catalog):
    lookup_dec = np.array([r["dec_deg"] for r in catalog])
    cases = []
    for field in evidence["fields"]:
        for source in field["sources"]:
            if not source["candidate"]:
                continue
            p = source["discovery"]["period_days"]
            nearby = []
            for j in np.flatnonzero(abs(lookup_dec-source["dec_deg"]) < RADIUS/3600):
                row = catalog[j]
                d = separation(source["ra_deg"], source["dec_deg"], row["ra_deg"], row["dec_deg"])
                if d > RADIUS:
                    continue
                test = [(abs(p/(tim*factor)-1), factor, tim)
                        for tim in row.get("period_days",[]) for factor in (0.5, 1., 2.)]
                order = min(test) if test else None
                nearby.append({**row, "separation_arcsec": d,
                               "best_fractional_period_error": order[0] if order else None,
                               "period_scale": order[1] if order else None,
                               "consistent": bool(order and order[0] <= PERIOD_TOLERANCE)})
            nearby.sort(key=lambda r:r["separation_arcsec"])
            cases.append({"field_id": field["field_id"],"source_id":source["source_id"],
                          "image_period_days":p,"counterparts":nearby})
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--evidence",type=Path,default=Path("results/ffi-lmc-grid.json"))
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    records, receipts = tables(args.cache)
    raw = args.evidence.read_bytes()
    results = {"schema_version":1,"input_sha256":sha256(raw).hexdigest(), "catalogs":receipts,
               "radius_arcsec":RADIUS,"period_tolerance_fraction":PERIOD_TOLERANCE,
               "period_scales":[0.5,1,2],"identifier_rows":len(records),
               "cases":match(json.loads(raw),records)}
    atomic_bytes(args.output,(json.dumps(results,indent=2,allow_nan=False)+"\n").encode())


if __name__ == "__main__":
    main()
