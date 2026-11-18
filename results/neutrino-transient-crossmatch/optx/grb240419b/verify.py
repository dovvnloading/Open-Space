"""Recompute spherical offsets, delay, and the observed XRT time span."""
import csv
from datetime import datetime
from math import asin,cos,degrees,radians,sin,sqrt
from pathlib import Path

d=Path(__file__).resolve().parent
def angle(a,b):
    x,y=map(radians,a)
    z,t=map(radians,b)
    return degrees(2*asin(sqrt(sin((t-y)/2)**2+cos(y)*cos(t)*sin((z-x)/2)**2)))
nu=(327.08,3.06)
xrt=(325.89670,4.23468)
for name,point in [("Swift GRB table BAT",(325.853,4.299)),("BAT mission reanalysis",(325.8913,4.232397)),
                   ("final enhanced XRT",xrt),("Nordic Optical Telescope",(325.8970375,4.2348056))]:
    print(name,"offset to neutrino [deg]",round(angle(point,nu),5),"offset to XRT [arcmin]",round(angle(point,xrt)*60,4))
lag=(datetime.fromisoformat("2024-04-24T01:49:26Z")-datetime.fromisoformat("2024-04-19T12:22:52.33Z")).total_seconds()/86400
print("lag [d]",round(lag,7))
data=[tuple(map(float,l.strip().split())) for l in (d/"xrt_curve.qdp").read_text().splitlines() if l and l[0].isdigit()]
print("XRT bins",len(data),"last bin center [d]",max(x[0] for x in data)/86400,"last right edge [d]",max(x[0]+x[1] for x in data)/86400)
