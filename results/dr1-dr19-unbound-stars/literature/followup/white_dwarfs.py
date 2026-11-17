"""Cross-check the DESI high-speed selection with the DR1 white-dwarf specialist fits."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
from urllib.request import urlopen

SOURCE = "https://zenodo.org/api/records/23453342/files/DESI_DR1_catalogue_1.0.csv.gz/content"
HERE = Path(__file__).resolve().parent
KEEP = ["WDJname","specType","specType_confidence","model","model_confidence","RV","RV_error",
        "nExposures","primary","logP_Huber","distance","distance_err","mass","fluxContamination_DESI"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", type=Path)
    parser.add_argument("--catalogue", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE/"white_dwarf_matches.csv")
    args = parser.parse_args()
    if not args.catalogue.exists():
        args.catalogue.write_bytes(urlopen(SOURCE,timeout=180).read())
    if hashlib.sha256(args.catalogue.read_bytes()).hexdigest() != "b86d315dcee214d0fd49ecf1012f46e52306144a64da8baf079390b87403eaa5":
        raise ValueError("White-dwarf catalogue changed")
    with args.selection.open() as f:
        sources = {r["SOURCE_ID"]: r for r in csv.DictReader(f)}
    with gzip.open(args.catalogue,"rt") as f, args.output.open("w",newline="") as out:
        w=csv.DictWriter(out,fieldnames=["source_id","nearby"]+KEEP,lineterminator="\n")
        w.writeheader()
        for r in csv.DictReader(f):
            sid=r["designation"].split()[-1]
            if sid in sources:
                w.writerow({"source_id":sid,"nearby":sources[sid]["nearby"],**{k:r[k] for k in KEEP}})


if __name__ == "__main__":
    main()
