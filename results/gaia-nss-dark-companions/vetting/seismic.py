"""Replace the photometric primary mass by a published seismic mass, holding Gaia's mass function fixed."""
import argparse
import csv
import hashlib
import io
import json
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from scipy.optimize import brentq

SOURCES = {
    "apokasc4.txt": "https://content.cld.iop.org/journals/0067-0049/276/2/69/revision1/apjsad9feft4_mrt.txt",
    "kepler.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/696/A243/tablea1&-out=KIC,GaiaDR3&-out.max=unlimited",
    "apok2.csv": "https://raw.githubusercontent.com/Jesstella/APO-K2/c9fafef896cf82c352aba05aa627aefc6000d12d/apo_k2_public_catalog.csv",
}


def secondary(f, primary):
    return brentq(lambda m: m**3/(m+primary)**2-f, 0.0, 1000.0)


def main(cache, out):
    cache.mkdir(exist_ok=True, parents=True)
    root = Path(__file__).resolve().parents[1]
    inputs = {}
    for filename, url in SOURCES.items():
        path = cache / filename
        if not path.is_file():
            with urllib.request.urlopen(url, timeout=90) as response:
                path.write_bytes(response.read())
        inputs[filename] = dict(url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    table = {}
    for screen in ("extension/mass_screen.csv", "gaia/mass_screen.csv"):
        path = root/screen
        inputs[screen] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        with path.open() as stream:
            for r in csv.DictReader(stream):
                # The native DR3 binary-primary screen takes priority on overlap.
                key=r["source_id"]
                if key not in table or float(table[key]["m2_screen_msun"])<float(r["m2_screen_msun"]) and table[key]["screen"]==screen or table[key]["screen"]!=screen:
                    table[key] = dict(r,screen=screen)
    source_by_kic = defaultdict(set)
    kics_by_source = defaultdict(set)
    for line in (cache/"kepler.tsv").read_text().splitlines():
        parts = line.strip().split()
        if len(parts)==2 and all(p.isdecimal() for p in parts):
            kic, source = parts
            source_by_kic[kic].add(source)
            kics_by_source[source].add(kic)
    catalog_rows = []
    counts=Counter()
    for line in (cache/"apokasc4.txt").read_text().splitlines():
        kic = line[:8].strip()
        if not kic.isdecimal() or len(line)<148:continue
        counts["APOKASC3_table_rows"]+=1
        if len(source_by_kic[kic])!=1:continue
        sid=next(iter(source_by_kic[kic]))
        if sid not in table or len(kics_by_source[sid])!=1:continue
        counts["APOKASC3_screen_matches"]+=1
        label=line[22:30].strip();mass=float(line[126:136]);error=float(line[137:147])
        if label not in ("Gold","Silver") or not (mass>0 and error>0):continue
        catalog_rows.append((sid,kic,"APOKASC-3",label,mass,error))
    with (cache/"apok2.csv").open() as f:
        for r in csv.DictReader(f):
            counts["APOK2_table_rows"]+=1
            sid=r["gaia_source"]
            if sid not in table:continue
            counts["APOK2_screen_matches"]+=1
            mass,err=float(r["mass"]),float(r["mass_err"])
            if mass>0 and err>0:catalog_rows.append((sid,r["epic"],"APO-K2",r["ev"],mass,err))
    rows=[]
    for sid,identifier,catalog,flag,m1,error in catalog_rows:
        original=table[sid];f=float(original["f_msun"])
        rows.append(dict(source_id=sid,catalog=catalog,catalog_id=identifier,category=flag,
                         screen=original["screen"],nss_solution_type=original["nss_solution_type"],mass_basis=original['mass_basis'],
                         f_msun=f,screen_m1=original["m1"],screen_m2=original['m2_screen_msun'],
                         seismic_m1=m1,seismic_m1_err=error,seismic_m2=secondary(f,m1),
                         m2_at_m1_minus_error=secondary(f,max(0.001,m1-error)),m2_at_m1_plus_error=secondary(f,m1+error)))
    rows.sort(key=lambda r:(int(r["source_id"]),r["catalog"]))
    assert len({r['source_id'] for r in rows})==len(rows), "Duplicate seismic source requires explicit adjudication"
    with (out/"seismic_masses.csv").open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    counts.update(usable=len(rows),below_1p4=sum(r["seismic_m2"]<1.4 for r in rows),
                  above_3=sum(r["seismic_m2"]>3 for r in rows),selected_sources=len(table))
    (out/"seismic_depth.json").write_text(json.dumps(dict(counts=dict(counts),inputs=inputs),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("cache",type=Path);p.add_argument("out",type=Path)
    args=p.parse_args();main(args.cache,args.out)
