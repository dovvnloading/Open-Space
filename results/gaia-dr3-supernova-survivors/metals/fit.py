"""Measure DESI DR1 coadds of J1240+6710 and LP 93-21 reusing the J1109 fitter."""
import argparse
import csv
import hashlib
import io
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"spectra"))
from fit_control import C, data_window, ew_window, velocity_fit

OBJECTS = {"LP 93-21":"WDJ104559.42+590451.32","J1240+6710":"WDJ124043.46+671035.93"}
HASHES = {'WDJ124043.46+671035.93': '234869b34e7c4b794fe317f75a6c37793008893bff4c2e44c6ae89633236fcdb', 'WDJ104559.42+590451.32': 'e137e16bafab2d9daf959d033a1da2906483ce9ee3720749dcc4bf7dbb75c360'}


def carbon_shift(spectra, rest, half):
    w, f, weights = data_window(spectra,rest-half,rest+half)
    x = w-rest
    def residual(q):
        level,slope,depth,cen,width = q
        return (level+slope*x-depth*np.exp(-0.5*((x-cen)/width)**2)-f)*weights
    bounds = ([1,-10,0,-12,.45],[200,10,200,18,15])
    fits = [least_squares(residual,[np.median(f),0,8,cen,4],bounds=bounds,x_scale="jac")
            for cen in [-8,-2,1,4,7,12]]
    fit = min(fits,key=lambda r:r.cost)
    cov = np.linalg.pinv(fit.jac.T@fit.jac)*max(1,2*fit.cost/(len(w)-5))
    ratio = 1.+fit.x[3]/rest
    v = C*(ratio**2-1.)/(ratio**2+1.)
    e = 4*C*ratio/(rest*(ratio**2+1.)**2)*np.sqrt(cov[3,3])
    return dict(shift_kms=round(v,2),formal_sigma_kms=round(e,2),width_A=round(fit.x[4],2))


def write_csv(name, rows):
    with (HERE/name).open("w",newline="") as out:
        writer = csv.DictWriter(out,fieldnames=rows[0],lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(cache):
    cache.mkdir(parents=True,exist_ok=True)
    receipts, shifts, features = [], [], []
    for name, star in OBJECTS.items():
        url = f"https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/{star[2:7]}/{star}/{star}_DESI_DR1_spectra.tar.gz"
        path = cache/(star+".tar.gz")
        if not path.exists():
            with urllib.request.urlopen(url,timeout=90) as response:
                path.write_bytes(response.read())
        raw = path.read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        if sha != HASHES[star]:
            raise ValueError("Archive hash mismatch for "+star)
        arms, members = [], []
        with tarfile.open(fileobj=io.BytesIO(raw)) as tar:
            for m in tar.getmembers():
                if m.isfile() and m.name.endswith(".dat"):
                    b = tar.extractfile(m).read()
                    a = np.loadtxt(io.BytesIO(b))
                    if "_coadd-" in m.name:
                        arms.append(a)
                    members.append(dict(file=m.name,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),npix=len(a)))
        assert len(arms)==3 and len(members)==6
        receipts.append(dict(name=name,url=url,sha256=sha,bytes=len(raw),members=members))
        if name=="LP 93-21":
            # Laboratory vacuum wavelengths as tabulated by Kawka et al. (2020).
            for rest in [4933.40,5053.56,5381.83]:
                for half in [13,25,40]:
                    shifts.append(dict(name=name,region=f"C I vacuum {rest}",half_width_A=half,
                                       **carbon_shift(arms,rest,half)))
        else:
            # Mg II is an unresolved 4481 multiplet; 6158 and 7773 are O blends.
            for label, air in [("Mg II 4481",[4481.20]),("Si II 6347/6371",[6347.11,6371.37]),
                               ("O I 6158",[6158.19]),("O I 7773 blend",[7773.44])]:
                result = velocity_fit(arms,air)
                shifts.append(dict(name=name,region=label,half_width_A=13,shift_kms=round(result["rv"],2),
                                   formal_sigma_kms=round(result["formal_sigma"],2),width_A=""))
            # Externally fixed literature apparent photospheric shift; do not use
            # gravity-corrected systemic velocity to place absorption windows.
            for label, air in [("Mg II 4481",[4481.2]),("O I 6158",[6158.19]),("O I 7774",[7771.94,7774.17,7775.39]),
                               ("Ne I 6402",[6402.25]),("Si II 6347",[6347.11]),("Si II 6371",[6371.37]),
                               ("H alpha",[6562.80]),("H beta",[4861.33]),("He I",[5875.62])]:
                ew = ew_window(arms,air,-158.)
                features.append(dict(name=name,region=label,pseudo_EW_A=round(ew[0],3) if ew else "",
                                     formal_sigma_A=round(ew[1],3) if ew else "",apparent_shift_kms=-158.))
    (HERE/"spectral_receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")
    write_csv("spectral_shifts.csv",shifts)
    write_csv("feature_indices.csv",features)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("cache",type=Path)
    run(parser.parse_args().cache)
