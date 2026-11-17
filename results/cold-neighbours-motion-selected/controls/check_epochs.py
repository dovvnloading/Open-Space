"""Reassess control identities at the WISE epoch using independent census astrometry."""
import csv
import math
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"pilot"))
from select_catwise import select,load_corrections
EPOCH=2015.405

def position(a,d,mx,my,dt):
    return (a+dt*mx/(3600*math.cos(math.radians(d))))%360,d+dt*my/3600
def distance(a,d,b,e):
    return 3600*math.hypot(((a-b+180)%360-180)*math.cos(math.radians(d)),d-e)

def main():
    inputs=ROOT.parent/"catalogues"/"inputs"
    census=[]
    for r in csv.DictReader((inputs/"20pc").open()):
        try:
            mx=float(r["PMRA (mas/yr)"])/1000;my=float(r["PMDec (mas/yr)"])/1000
            a,d=position(float(r["RA (deg)"]),float(r["Dec (deg)"]),mx,my,EPOCH-float(r["Epoch"]))
            census.append((r["DefaultName"],a,d,mx,my,float(r["Parallax (mas)"])))
        except ValueError:continue
    ofs=load_corrections();outputs=[]
    for r in csv.DictReader((ROOT/"recovery.csv").open()):
        a,d,mx,my,plx=map(float,[r[k] for k in ["ra2015","dec2015","pmra","pmdec","parallax_mas"]])
        # Avoid unrelated companions or crossings. Prefer the exact designation when present.
        matches=[s for s in census if distance(a,d,s[1],s[2])<=12 and math.hypot(mx-s[3],my-s[4])<=max(.2,.3*math.hypot(mx,my)) and abs(plx-s[5])<=max(20,.2*plx)]
        exact=[s for s in matches if s[0].strip()==r["control"].strip()]
        ref=(exact or matches)
        anchor="UltracoolSheet J2000 propagated"; shift=0
        if len(ref)==1:
            n,b,e,mx,my,p=ref[0];anchor="20pc: "+n;shift=distance(a,d,b,e);a,d=b,e
        data=list(csv.DictReader((ROOT/"inputs"/r["file"]).open()));counts={2.5:[],5.:[],10.:[]}
        for z in data:
            row=select(z,ofs)
            if not row:continue
            if math.hypot(row["pmra"]-mx,row["pmdec"]-my)>max(.2,.3*math.hypot(mx,my)):continue
            sep=distance(a,d,row["ra2015_deg"],row["dec2015_deg"])
            for limit in counts:
                if sep<=limit:counts[limit].append((sep,row))
        out=dict(control=r["control"],parallax_mas=plx,anchor=anchor,shift_arcsec=round(shift,3),ra2015=a,dec2015=d)
        for limit,ls in counts.items():
            tag={2.5:"2p5",5.:"5",10.:"10"}[limit];ls.sort(key=lambda v:v[0])
            out.update({f"rows_{tag}":len(ls),f"bright_rows_{tag}":sum(x[1]["bright_tier"] for x in ls),
                        f"closest_{tag}":ls[0][1]["source_name"] if ls else ""})
        outputs.append(out)
    with (ROOT/"epoch_matched.csv").open("w") as f:
        writer=csv.DictWriter(f,fieldnames=outputs[0]);writer.writeheader();writer.writerows(outputs)
    print(len(outputs),sum(x['anchor'].startswith('20pc') for x in outputs))
    for k in ['rows_2p5','bright_rows_2p5','rows_5','bright_rows_5','rows_10']:print(k,sum(x[k]>0 for x in outputs))
    for x in outputs:
        if x['parallax_mas']>=200:print(x)
if __name__=="__main__":main()
