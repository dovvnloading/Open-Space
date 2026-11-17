"""Retrieve RVSpecFit atmosphere parameters for temperature-stratified false-positive rates."""
import argparse
import csv
from pathlib import Path
from catalog_audit import DESI, freeze


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("roster",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    with a.roster.open() as f:
        pixels=sorted({int(r["healpix"]) for r in csv.DictReader(f)})
    columns="targetid,healpix,teff,logg,feh,sn_b,sn_r,sn_z,rvs_warn,success"
    sql="SELECT "+columns+" FROM desi_dr1.mws WHERE survey='main' AND program='backup' AND healpix IN ("
    a.out.mkdir(exist_ok=True,parents=True)
    for n in range(0,len(pixels),12):
        groups=",".join(map(str,pixels[n:n+12]))
        freeze(a.out,f"parameters_{n//12:02d}",DESI,sql+groups+")",columns.split(","))
