"""Retrieve merged-source Euclid PHZ for AstroVink positions with a subarcsecond match."""
import argparse
import csv
import json
from pathlib import Path
from fetch_phz import ajax, FIELDS


def main(cache, out):
    with (out/"astrovink_mer_links.csv").open() as stream:
        ids=sorted({r["closest_mer_id_within_0p5arcsec"] for r in csv.DictReader(stream)
                    if r["closest_mer_id_within_0p5arcsec"]},key=int)
    path=cache/"phz_batches"
    path.mkdir(exist_ok=True,parents=True)
    rows,manifest=ajax(ids,path)
    with (out/"astrovink_phz.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS,lineterminator="\n");writer.writeheader()
        writer.writerows(sorted(rows,key=lambda r:int(r["object_id"])))
    (out/"astrovink_phz_query.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(len(rows),"rows for",len(ids),"centres")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir",type=Path)
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args(); main(a.input_dir,a.output)
