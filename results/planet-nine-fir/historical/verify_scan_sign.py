"""Verify the orientation of the SCANPI arcminute axis using a second sky origin."""
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
from urllib.request import urlopen

import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent / "scanpi"


def read(raw):
    lines = raw.decode().splitlines()
    curve = np.array([[float(s) for s in line.split()[:3]] for line in lines if line.startswith(" ")])
    info = dict(re.findall(r"^\\(\w+)\s*=\s*([^\n]+)",raw.decode(),re.M))
    return curve[:,0],curve[:,2],info


def main():
    examples=json.loads((ROOT/'sign_inputs.json').read_text())
    rows=[]
    with tarfile.open(ROOT/'profiles.tar.gz') as fsc, tarfile.open(ROOT/'galaxy_profiles.tar.gz') as galaxy:
        for item in examples:
            member=item['name']
            x,f,orig=read(fsc.extractfile(member).read())
            u,g,gal=read(galaxy.extractfile(member).read())
            sample=abs(x)<30
            def loss(offset):
                diff=f[sample]-np.interp(x[sample]-offset,u,g)
                return np.mean((diff-np.median(diff))**2)
            best=minimize_scalar(loss,bounds=(-.75,.75),method='bounded').x
            pa=np.deg2rad(float(orig['pa']))
            # Tangent-plane east/north differences, in arcmin, accurate to <0.01 arcsec here.
            east=(35.7544007-35.74075)*60*np.cos(np.deg2rad(-48.5125))
            north=(-48.5112676+48.5125)*60
            along=east*np.sin(pa)+north*np.cos(pa)
            cross=east*np.cos(pa)-north*np.sin(pa)
            rows.append(dict(member=member,pa_deg=float(orig['pa']),rev=int(orig['rev']),
                  predicted_galaxy_along_arcmin=along,measured_galaxy_along_arcmin=best,
                  predicted_galaxy_cross_arcmin=cross,
                  galaxy_minus_fsc_cross_scan_arcmin=float(gal['cross_scan'])-float(orig['cross_scan'])))
    with (ROOT/'sign_check.csv').open('w',newline='') as file:
        writer=csv.DictWriter(file,rows[0]);writer.writeheader();writer.writerows(rows)


if __name__=='__main__': main()
