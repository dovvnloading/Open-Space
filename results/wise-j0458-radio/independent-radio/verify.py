"""Fixed-beam photometry on independent low-frequency survey products."""
import csv
import hashlib
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
STAR = (74.7268, 64.5825)
PM = (218., 300.)
EPOCH = 60859.


def position(mjd):
    years = (mjd - EPOCH) / 365.25
    return (STAR[0] + PM[0]*years/3.6e6/np.cos(np.deg2rad(STAR[1])),
            STAR[1] + PM[1]*years/3.6e6)


def fit_at(img, x, y, jac, axes):
    ya, xa = np.indices(img.shape)
    en = jac @ np.stack([xa-x, ya-y]).reshape(2,-1)
    east, north = en.reshape(2,*img.shape)
    g = np.exp(-4*np.log(2)*((east/axes[0])**2+(north/axes[1])**2))
    take = (g>1e-5) & np.isfinite(img)
    basis = np.column_stack([g[take],np.ones(take.sum()),east[take],north[take]])
    return np.linalg.lstsq(basis,img[take],rcond=None)[0][0]


def main():
    with (ROOT/"sources.csv").open() as f:
        sources = list(csv.DictReader(f))
    results = []
    # Evaluate all surveys on identical sightlines. Only the 2023 anchor
    # approximates a LOFAR observation date; NVSS/WENSS have much bigger beams.
    with (ROOT/"lotss-source-cone.csv").open() as f:
        cat = next(r for r in csv.DictReader(f) if r["source_name"]=="ILTJ045854.31+643455.8")
    coordinates = {"catalog":(float(cat["ra"]),float(cat["dec"])),
                   "star_2017":position(58037.288),
                   "star_2023":position(59970.715972),
                   "star_2025":position(60924.590542),
                   "star_1995":position(49800.)}
    for source in sources:
        name = source["local_product"]
        path = ROOT/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=source["sha256"]:
            raise ValueError(f"Input hash differs for {name}")
        with fits.open(path) as hdus:
            h = hdus[0].header
            d = hdus[0].data.squeeze().astype(float)*1e3
        w=WCS(h).celestial
        jac = w.pixel_scale_matrix*3600
        axes = (float(source["beam_east_arcsec"]),float(source["beam_north_arcsec"]))
        x0,y0=w.world_to_pixel_values(*coordinates["catalog"])
        step=int(np.ceil(1.5*axes[0]/np.sqrt(abs(np.linalg.det(jac)))))
        controls=[fit_at(d,x0+step*i,y0+step*j,jac,axes)
                  for i in range(-5,6) for j in range(-5,6) if 4<i*i+j*j<27]
        md=np.median(controls)
        spread=1.4826*np.median(np.abs(np.asarray(controls)-md))
        row=dict(product=name,n_controls=len(controls),control_median_mjy=md,
                 control_mad_scale_mjy=spread)
        for key,pos in coordinates.items():
            x,y=w.world_to_pixel_values(*pos)
            row[key+"_mjy"]=fit_at(d,x,y,jac,axes)
        if name.startswith("lotss"):
            # Independent circular-restoring-beam centroid on central patch.
            ygrid,xgrid=np.indices(d.shape)
            cut=np.hypot(xgrid-x0,ygrid-y0)<8
            e,n=(jac@np.array([xgrid[cut]-x0,ygrid[cut]-y0]))
            def res(p):
                amp,ce,cn,bg=p
                g=np.exp(-4*np.log(2)*((e-ce)**2+(n-cn)**2)/axes[0]**2)
                return amp*g+bg-d[cut]
            amp,ce,cn,bg=least_squares(res,[.9,0.,0.,0.]).x
            xo,yo=np.linalg.solve(jac,[ce,cn])
            ra,dec=w.pixel_to_world_values(x0+xo,y0+yo)
            row.update(free_peak_mjy=amp,fit_ra=float(ra),fit_dec=float(dec))
        results.append(row)
    with (ROOT/"independent-measurements.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(results)


if __name__ == "__main__":
    main()
