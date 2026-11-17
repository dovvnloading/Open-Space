"""Apply the same photometry/motion thresholds to blind fields and controls."""
import csv
import math
from collections import Counter
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
OFFSETS_URL = "https://irsa.ipac.caltech.edu/data/WISE/CatWISE/gator_docs/CatWISE2020_Table1_20201012.tbl"
FLOOR = 0.035  # arcsec/year in each axis, for ranking motion estimates


def number(row, key):
    return float(row[key]) if row[key].strip() else math.nan

def load_corrections():
    path = ROOT/"inputs"/"catwise_tile_offsets.tbl"
    if not path.exists():
        with urlopen(OFFSETS_URL,timeout=90) as f:
            path.write_bytes(f.read())
    offsets = {}
    for line in path.read_text().splitlines():
        if line.startswith(("|","\\")) or not line.strip(): continue
        tile, ra_true, _, dec, pmra, pmdec = line.split()
        offsets[tile] = tuple(map(float, (ra_true, dec, pmra, pmdec)))
    return offsets

def prefilter(r):
    """Same inclusion rule as the SQL, including W1 non-detections."""
    sn1 = number(r,"w1snr_pm")
    raw_mu = math.hypot(number(r,"pmra"),number(r,"pmdec"))
    return (number(r,"w2snr_pm") >= 8 and raw_mu >= .2
            and (not math.isfinite(sn1) or sn1 < 3
                 or number(r,"w1mpro_pm")-number(r,"w2mpro_pm") >= 1.5))

def select(r, offsets):
    return measure(r, offsets) if prefilter(r) else None

def measure(r, offsets):
    dpra, dpdec, dmra, dmdec = offsets[r["source_id"][:8]]
    ma = number(r,"pmra") + dmra
    md = number(r,"pmdec") + dmdec
    mu = math.hypot(ma,md)
    if not all(number(r,k) > 0 for k in ("sigpmra","sigpmdec")): return
    z = max(abs(ma)/math.hypot(number(r,"sigpmra"),FLOOR),
            abs(md)/math.hypot(number(r,"sigpmdec"),FLOOR))
    strong_w1 = number(r,"w1snr_pm") >= 3
    if strong_w1:
        color = number(r,"w1mpro_pm") - number(r,"w2mpro_pm")
    else:
        # One-sided 95% Gaussian bound per flux. A positive measured W1 flux
        # raises its upper bound; a negative measurement is replaced by zero.
        if not all(math.isfinite(number(r,k)) for k in
                   ("w1flux_pm","w1sigflux_pm","w2flux_pm","w2sigflux_pm")): return
        if min(number(r,"w1sigflux_pm"),number(r,"w2sigflux_pm")) <= 0: return
        w1hi=max(0,number(r,"w1flux_pm")) + 1.645*number(r,"w1sigflux_pm")
        w2lo=number(r,"w2flux_pm") - 1.645*number(r,"w2sigflux_pm")
        if w2lo <= 0: return
        color=2.5*math.log10(w2lo/w1hi)
    if not (color >= 2 and mu >= .3 and z >= 5): return
    # A motion component pinned to +/-100 arcsec/yr indicates a boundary fit.
    motion_boundary = max(abs(number(r,"pmra")),abs(number(r,"pmdec"))) >= 99.99
    w2cc = r["cc_flags"][1:2]
    bright_tier = (not motion_boundary and number(r,"w2mpro_pm") <= 16.5
                   and number(r,"w2rchi2_pm") <= 3 and r["w2ab_map"] == "0"
                   and w2cc in ("","0"))
    dec = number(r,"dec_pm")
    return dict(source_name=r["source_name"],source_id=r["source_id"],
        ra2015_deg=round(number(r,"ra_pm")+dpra/math.cos(math.radians(dec)),7),
        dec2015_deg=round(dec+dpdec,7),pmra=round(ma,6),pmdec=round(md,6),
        sigpmra=r["sigpmra"],sigpmdec=r["sigpmdec"],mu=round(mu,5),z_axis=round(z,3),
        w1w2=round(color,4),w1_mode="measured" if strong_w1 else "lower_bound",
        w2=r["w2mpro_pm"],sn2=r["w2snr_pm"],rchi2=r["w2rchi2_pm"],
        w2ab_map=r["w2ab_map"],w2cc=w2cc,pmcode=r["pmcode"],
        motion_boundary=int(motion_boundary),bright_tier=int(bright_tier))

def main():
    offsets = load_corrections()
    rows=[]
    for file in sorted((ROOT/"inputs").glob("*_rows.csv")):
        for r in csv.DictReader(file.open()):
            candidate=select(r,offsets)
            if candidate: rows.append(dict(field=file.stem.split("_")[0],**candidate))
    if not rows: print("zero rows"); return
    with (ROOT/"shortlist.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    print(Counter(r["field"] for r in rows))
    clean=[r for r in rows if float(r["rchi2"]) <= 3 and r["w2ab_map"]=="0" and r["w2cc"] in ("0","")]
    print("Diagnostic priority",len(clean))
    for r in sorted(clean,key=lambda r:-r["w1w2"])[:14]:
        print(r["field"],r["source_name"],r["w2"],r["w1w2"],r["mu"],r["rchi2"],r["w1_mode"])
if __name__=="__main__": main()
