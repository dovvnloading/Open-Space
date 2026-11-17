"""Measure the single-exposure blend against the track in fixed sky coordinates."""
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
RA,DEC,MA,MD=274.5676647347,66.42021838013,.055676706,-1.025546534
def main():
    groups=defaultdict(list)
    with (HERE/"l1b.csv").open() as f:
        for d in csv.DictReader(f):
            if float(d["w2snr"])<5 or int(d["qual_frame"])<5: continue
            x=(float(d["ra"])-RA)*math.cos(math.radians(DEC))*3600
            y=(float(d["dec"])-DEC)*3600
            if math.hypot(x,y)>13: continue
            mjd=float(d["mjd"]); dt=(mjd-57170)/365.25
            offset=math.hypot(x-MA*dt,y-MD*dt)
            groups[round((mjd-56647)/182.625)].append((mjd,x,y,offset,float(d["w2mpro"])))
    out=[]
    for season,rows in sorted(groups.items()):
        b=[s for s in rows if math.hypot(s[1]-.15,s[2]+2.95)<=2]
        out.append(dict(season=season,mjd=round(st.median(s[0] for s in rows),3),
          n_within_13=len(rows),n_within_1p5_of_track=sum(s[3]<=1.5 for s in rows),
          n_within_2_of_blend=len(b),blend_east=round(st.median(s[1] for s in b),3),
          blend_north=round(st.median(s[2] for s in b),3),blend_w2=round(st.median(s[4] for s in b),3)))
    with (HERE/"visits.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
if __name__=="__main__":main()
