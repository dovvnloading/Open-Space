"""Locate the 79-minute signal in the December 2023 MASTER J055845 field.

Run from the result directory after downloading the inputs described in README.
All pixel coordinates are zero based in the 7 by 7 TESSCut product.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import ndtr

FREQ = 18.179  # cycles/day, independent published Sector 73 frequency
BLOCK = 0.25
WINDOWS = {"pre_erupt": (3287, 3290.75), "early": (3292.5, 3299.5),
           "late": (3306.5, 3312.5)}
TARGET = "3458275544681382912"
EB = "3458275643465151744"


def load_cube(path):
    with fits.open(path) as hdus:
        data = hdus[1].data
        t, cube = data["TIME"].astype(float), data["FLUX"].astype(float)
        good = (data["QUALITY"] == 0) & np.isfinite(t) & np.isfinite(cube).all(axis=(1, 2))
        return t[good], cube[good], WCS(hdus[2].header)


def sincos(t, frequency):
    phase = 2 * np.pi * frequency * (t - 3293)
    return np.stack([np.sin(phase), np.cos(phase)], axis=1)


def block_fit(t, flux, start, frequency):
    use = (t >= start) & (t < start + BLOCK)
    if use.sum() < 60:
        return None
    u = (t[use] - start - BLOCK / 2) / BLOCK
    design = np.column_stack([np.ones_like(u), u, u**2, sincos(t[use], frequency)])
    coeff = np.linalg.lstsq(design, flux[use], rcond=None)[0]
    return coeff[-2:]


def fit_center(image):
    """Integrate an elliptical Gaussian over pixel footprints with a tilted background."""
    yy, xx = np.indices(image.shape)
    def residual(p):
        total, xc, yc, sx, sy, base, gx, gy = p
        px = ndtr((xx+.5-xc)/sx) - ndtr((xx-.5-xc)/sx)
        py = ndtr((yy+.5-yc)/sy) - ndtr((yy-.5-yc)/sy)
        return (total*px*py + base + gx*(xx-3) + gy*(yy-3) - image).ravel()
    pars = least_squares(residual, [max(np.ptp(image)*4,.01),3.5,3.2,.6,.6,np.median(image),0,0],
                         bounds=([0,1,1,.2,.2,-np.inf,-np.inf,-np.inf],
                                 [np.inf,5,5,2,2,np.inf,np.inf,np.inf]))
    return pars.x[:5]


def window(t, cube, lc, label, rng, block_rows):
    low, high = WINDOWS[label]
    starts = np.arange(low, high, BLOCK)
    data, scalar = [], []
    for s in starts:
        coeff = block_fit(t, cube.reshape(len(t),-1), s, FREQ)
        if coeff is None:
            continue
        c = block_fit(t, lc, s, FREQ)
        data.append(coeff.reshape(2,7,7))
        scalar.append(c)
        block_rows.append(dict(window=label,btjd=float(s+BLOCK/2),sin=float(c[0]),cos=float(c[1])))
    scalar, data = np.array(scalar), np.array(data)
    mean = scalar.mean(axis=0)
    # A single rotation per window preserves the relative phases of all pixels.
    phase = mean / np.linalg.norm(mean)
    signed = np.einsum("bijk,i->bjk", data, phase)
    image = signed.mean(axis=0)
    # Resampling the same block preserves contemporaneous correlations across 49 pixels.
    picks = rng.integers(len(signed), size=(300,len(signed)))
    centroids = np.array([fit_center(signed[g].mean(axis=0))[1:3] for g in picks])
    obs = np.linalg.norm(mean)
    side = []
    for _ in range(50):
        angles = rng.uniform(0,2*np.pi,size=(1000,len(scalar)))
        z = scalar[:,0]+1j*scalar[:,1]
        side.extend(abs((np.exp(1j*angles)*z).mean(axis=1)))
    fp = float((1+np.sum(np.asarray(side)>=obs))/(1+len(side)))
    return {"blocks":len(scalar), "semiamplitude_aperture_e_s":float(obs),
            "phase_cycles":float(np.arctan2(mean[1],mean[0])/(2*np.pi)),
            "block_phase_randomization_p":fp,
            "center_xy":fit_center(image)[1:3].tolist(),
            "block_bootstrap_xy_percentiles":np.percentile(centroids,[16,50,84],axis=0).tolist(),
            "signed_image_e_s":image.round(5).tolist()}


def epoch_period(t, lc, low, high):
    mask = (t>=low)&(t<high)
    tt,y=t[mask],lc[mask]
    residual=y-median_filter(y, size=73, mode="nearest")
    def power(f):
        a=sincos(tt,f)
        v=np.linalg.lstsq(a,residual,rcond=None)[0]
        return np.var(a@v)/np.var(residual)
    fgrid=np.arange(16,20,.015)
    best=fgrid[np.argmax([power(f) for f in fgrid])]
    f=minimize_scalar(lambda v:-power(v),bounds=(best-.02,best+.02),method="bounded").x
    return dict(n=len(tt),frequency_cpd=f,period_min=1440/f,fraction_variance=power(f))


def main():
    p=argparse.ArgumentParser()
    p.add_argument("fits");p.add_argument("gaia")
    p.add_argument("--output",default="measurement.json")
    p.add_argument("--blocks",default="blocks.csv")
    args=p.parse_args()
    t,cube,w=load_cube(args.fits)
    g=pd.read_csv(args.gaia,dtype={"source_id":str})
    # First-order Gaia propagation, with pmra = d(alpha)/dt * cos(delta).
    dt=2023.95-2016
    ra=g.ra+dt*g.pmra.fillna(0)/3.6e6/np.cos(np.deg2rad(g.dec))
    dec=g.dec+dt*g.pmdec.fillna(0)/3.6e6
    g["x"],g["y"]=w.world_to_pixel_values(ra,dec)
    target=g[g.source_id==TARGET].iloc[0]
    yy,xx=np.indices((7,7)); r2=(xx-target.x)**2+(yy-target.y)**2
    # Use the catalog position to define the aperture before looking at its clock.
    weights=np.exp(-r2/(2*.80**2))*(r2<4)
    bg=np.median(cube[:,r2>=9],axis=1)
    lc=np.einsum("tij,ij->t",cube-bg[:,None,None],weights)
    rng=np.random.default_rng(47173)
    out={"quality_zero_finite_cadences":len(t),"input_frequency_cpd":FREQ,"windows":{},
         "nearby_gaia":[r for r in g.loc[:,
                    ["source_id","x","y","phot_g_mean_mag"]].to_dict("records")
                    if np.hypot(r["x"]-target.x,r["y"]-target.y)<1.4]}
    rows=[]
    for label,(lo,hi) in WINDOWS.items():
        out["windows"][label]=window(t,cube,lc,label,rng,rows)
        out["windows"][label]["frequency_fit"]=epoch_period(t,lc,lo,hi)
    # Localize the rapid eruption rise against a precursor window on the same TESS orbit.
    def median_image(lo,hi):
        return np.median(cube[(t>=lo)&(t<hi)],axis=0)
    out["rise_difference_centers"]=[{"btjd":v,"xy":fit_center(median_image(*v[2:])-median_image(*v[:2]))[1:3].tolist()}
        for v in [[3290.6,3290.9,3291.3,3291.6],[3290.5,3290.9,3291.5,3291.9],
                  [3290,3290.4,3291.4,3291.8]]]
    pd.DataFrame(rows).to_csv(args.blocks,index=False,float_format="%.8f")
    Path(args.output).write_text(json.dumps(out,indent=2,allow_nan=False,default=float)+"\n")


if __name__=="__main__":
    main()
