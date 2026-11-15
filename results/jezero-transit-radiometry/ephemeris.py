"""Recalculate geometric screen with NAIF/JPL geometry, independently of the radiometry."""
import argparse
from pathlib import Path
import re
import urllib.request

import numpy as np
import pandas as pd
import spiceypy as sp
from spiceypy.utils.support_types import SPICEDOUBLE_CELL

ROOT = "https://naif.jpl.nasa.gov/pub/naif/"
SOL_BOUNDARIES = [0,89,179,299,419,539,659,779,899,1019,1139,1259,1379,1499,1619,1739,1859]


def load_kernels(cache):
    cache.mkdir(parents=True,exist_ok=True)
    files = ["generic_kernels/lsk/naif0012.tls", "MARS2020/kernels/pck/pck00010.tpc",
             "MARS2020/kernels/sclk/M2020_168_SCLKSCET.00009.tsc", "MARS2020/kernels/spk/de438s.bsp",
             "generic_kernels/spk/satellites/mar099s.bsp"]
    for start,end in zip(SOL_BOUNDARIES,SOL_BOUNDARIES[1:]):
        files.append(f"MARS2020/kernels/spk/m2020_surf_rover_loc_{start:04d}_{end:04d}_v1.bsp")
    for key in files:
        path=cache / Path(key).name
        if not path.exists():
            urllib.request.urlretrieve(ROOT+key,path)
        sp.furnsh(str(path))


def conjunctions():
    interval=SPICEDOUBLE_CELL(4)
    sp.wninsd(sp.str2et("2021-02-19"),sp.str2et("2026-05-12"),interval)
    # The search is over angular-separation minima, independent of how short optical events are.
    found=sp.gfsep("SUN","POINT","J2000","PHOBOS","POINT","J2000","LT+S","-168",
                   "LOCMIN",0,0,1500,20000,interval,SPICEDOUBLE_CELL(30000))
    anchor=sp.scs2e(-168,"1/780969452:00000") # archived LMST 01285M00:00:01.244
    rows=[]
    for idx in range(sp.wncard(found)):
        t=sp.wnfetd(found,idx)[0]
        sun=sp.spkpos("SUN",t,"IAU_MARS","LT+S","-168")[0]
        moon=sp.spkpos("PHOBOS",t,"IAU_MARS","LT+S","-168")[0]
        sep=sp.vsep(sun,moon)
        solar_radius=np.arcsin(695700/sp.vnorm(sun))
        bound=np.arcsin(13.4/sp.vnorm(moon))
        if sep>5*solar_radius:
            continue
        outward=sp.spkpos("-168",t,"IAU_MARS","NONE","499")[0]
        altitude=90-np.degrees(sp.vsep(sun,outward))
        if altitude<=5:
            continue
        clock=int(re.split(r"[-.:]",sp.sce2s(-168,t).split("/")[1])[0])
        sol=int(np.floor(1285+(t-anchor+1.244*1.02749)/88775.244))
        clearance=sep-solar_radius-bound
        rows.append(dict(sol=sol,moon="PHOBOS",utc=sp.et2utc(t,"ISOC",3),sclk=clock,alt=round(altitude,2),
                         separation_solar_radii=round(sep/solar_radius,4),
                         moon_solar_radii_max=round(bound/solar_radius,4),
                         geometry="off_disk" if clearance>=0 else "possible_disk_overlap"))
    return pd.DataFrame(rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,default=Path("kernels"))
    p.add_argument("--output",type=Path,default=Path("geometry.csv"))
    args=p.parse_args()
    load_kernels(args.cache)
    conjunctions().to_csv(args.output,index=False)
