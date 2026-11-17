"""Recover S/2025 U1's sky position with URA184 and check against Cordelia."""
import argparse
import csv
from pathlib import Path
import numpy as np
from astropy.time import Time
from geometry import sky_xy
from moon_spk import offsets


def run(directory, kernel):
    directory = Path(directory)
    with (directory/"selection.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["target"]=="URANUS"]
    with (directory/"ephemerides.csv").open() as f:
        eph = {(r["rootname"],int(r["body"])):r for r in csv.DictReader(f)}
    dates = np.array([(float(r["mjd_start"])+float(r["mjd_end"]))/2 for r in rows])
    ids = [r["rootname"] for r in rows]
    body = np.array([[float(eph[k,799][s]) for s in ("ra","dec","phase_deg","r_au","delta_au")] for k in ids])
    epoch = Time(dates,format="mjd",scale="utc").tdb.jd
    corrected = {}
    for target in (706,75052):
        vec = offsets(kernel,epoch,body,target)
        ra,dec = np.rad2deg(np.arctan2(vec[:,1],vec[:,0]))%360, np.rad2deg(np.arctan2(vec[:,2],np.hypot(vec[:,0],vec[:,1])))
        corrected[target] = ra,dec
    expected = [[float(eph[k,706][s]) for k in ids] for s in ("ra","dec")]
    residual = sky_xy(*corrected[706],*expected,1.,(0.,0.))
    err = np.hypot(*residual.T)
    if np.max(err) > .01:
        raise RuntimeError(f"Cordelia cross-check failed: max {np.max(err):.4f} arcsec")
    with (directory/"s2025u1_spk.csv").open("w",newline="") as f:
        out = csv.writer(f)
        out.writerow(["rootname","mjd_mid_utc","s2025u1_ra_deg","s2025u1_dec_deg","cordelia_spk_minus_horizons_mas"])
        for k,t,ra,dec,v in zip(ids,dates,*corrected[75052],err):
            out.writerow([k,f"{t:.10f}",f"{ra:.9f}",f"{dec:.9f}",f"{1000*v:.4f}"])
    print(f"{len(ids)} positions. Cordelia median/max residual {np.median(err)*1000:.3f}/{np.max(err)*1000:.3f} mas.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory")
    p.add_argument("kernel")
    a=p.parse_args()
    run(a.directory,a.kernel)
