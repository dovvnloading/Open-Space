"""Get Q1 MER central aperture diagnostics for candidate object IDs and AstroVink positions."""
import argparse
import csv
import hashlib
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from fetch_phz import ENDPOINT
from screen import output_csv, separation

BANDS = ["vis","y","j","h"]
FIELDS = ["object_id","right_ascension","declination","det_quality_flag","deblended_flag"]
for band in BANDS:
    FIELDS.extend([f"flux_{band}_1fwhm_aper",f"fluxerr_{band}_1fwhm_aper",f"flag_{band}"])


def query(where, cache):
    adql = "SELECT "+",".join(FIELDS)+" FROM q1.mer_catalogue WHERE "+where+" ORDER BY object_id"
    path = cache / (hashlib.sha256(adql.encode()).hexdigest()[:20]+".csv")
    if not path.exists():
        info = urlencode({"LANG":"ADQL","REQUEST":"doQuery","FORMAT":"csv","QUERY":adql})
        for attempt in range(4):
            try:
                with urlopen(ENDPOINT+"?"+info, timeout=100) as stream:
                    content = stream.read()
                break
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
    else:
        content = path.read_bytes()
    reader = csv.DictReader(io.StringIO(content.decode()))
    if set(reader.fieldnames or []) != set(FIELDS):
        raise ValueError("Unexpected MER columns: "+content[:250].decode())
    rows = list(reader)
    path.write_bytes(content)
    return rows, {"query": adql, "sha256":hashlib.sha256(content).hexdigest(), "rows":len(rows)}


def main(directory, out):
    def table(name):
        with (out/name).open() as stream:
            return list(csv.DictReader(stream))
    ids = {r["object_id"] for r in [*table("slde_screen.csv"),*table("agilelens_screen.csv")] if r["object_id"]}
    astro = table("astrovink_screen.csv")
    # The published tables sometimes center different catalogue entries on the same sky system.
    ids.update(r["archive_nearest_ab_id"].split("_")[-1] for r in astro if float(r["archive_nearest_ab_arcsec"])<10)
    ids = sorted(ids,key=int)
    conditions = ["object_id IN ("+",".join(ids[i:i+55])+")" for i in range(0,len(ids),55)]
    # Direct sub-arcsec search at paper coordinates also covers entries omitted from the archive.
    unmatched = [r for r in astro if not r["zenodo_id"]]
    terms = [f"1=CONTAINS(POINT('ICRS',right_ascension,declination),CIRCLE('ICRS',{r['ra_deg']},{r['dec_deg']},0.00014))" for r in unmatched]
    conditions += [term for term in terms]
    cache = directory/"mer_batches"
    cache.mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(lambda clause: query(clause,cache),conditions))
    combined = {r["object_id"]:r for block,_ in results for r in block}
    records = []
    for key in sorted(combined, key=int):
        r=combined[key]
        outrow={k:r[k] for k in FIELDS}
        for b in BANDS:
            f,e=r["flux_"+b+"_1fwhm_aper"],r["fluxerr_"+b+"_1fwhm_aper"]
            outrow[b+"_snr"]=float(f)/float(e) if f and e and float(e)>0 else ""
        records.append(outrow)
    output_csv(out,"mer_apertures.csv",records)
    links=[]
    for src in astro:
        possible=sorted((separation(float(src["ra_deg"]), float(src["dec_deg"]),float(r["right_ascension"]),
                                         float(r["declination"])),r["object_id"]) for r in records)
        d,key=possible[0]
        links.append({"name":src["name"],"grade":src["grade"],"closest_mer_id_within_0p5arcsec":key if d<0.5 else "",
                      "separation_arcsec":d if d<0.5 else ""})
    output_csv(out,"astrovink_mer_links.csv",links)
    (out/"mer_queries.json").write_text(json.dumps([manifest for _,manifest in results],indent=2)+"\n")
    print(f"{len(records)} distinct central MER entries, {sum(bool(r['closest_mer_id_within_0p5arcsec']) for r in links)}/36 AstroVink paper centers")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir",type=Path)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args()
    main(args.input_dir,args.output)
