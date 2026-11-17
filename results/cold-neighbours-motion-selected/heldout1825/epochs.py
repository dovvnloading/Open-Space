"""Native W2 same-season displacement test for the CatWISE J182532.09+671906.2 track."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import csv
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates

BASE = Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent/"control0855"))
from native import cutout
RA, DEC, PMX, PMY = 276.3837701, 67.3183865, -2.142453, 1.619898
CD = np.cos(np.deg2rad(DEC))
SCALE=1.375
# Choose two same-season visits, ten years apart, and the central 64 quality-10 exposures in each.
WINDOWS=[(56658.,56678.),(60288.,60310.)]


def build_manifest(ra, dec, windows, count, quality):
    rows=[]
    for group,(start,end) in enumerate(windows):
        url="https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?"+urlencode(dict(
            POS=f"{ra},{dec}", SIZE="0.001", WHERE=f"band=2 and mjd_obs>{start} and mjd_obs<{end} and qual_frame>={quality}",
            COLUMNS="scan_id,frame_num,mjd_obs,crval1,crval2,moon_sep,saa_sep"))
        archive=Table.read(urlopen(url,timeout=120).read().decode().splitlines(),format="ascii.ipac").to_pandas()
        archive["r"]=np.hypot((archive.crval1-ra)*np.cos(np.deg2rad(dec)),archive.crval2-dec)
        for f in archive.nsmallest(count,"r").sort_values("mjd_obs").itertuples():
            sid,frame=str(f.scan_id),f"{f.frame_num:03d}"
            name=f"{sid}{frame}-w2-int-1b.fits"
            img=("https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/"
                 f"p1bm_frm/{sid[-2:]}/{sid}/{frame}/{name}")
            rows.append(dict(epoch=group,mjd=f.mjd_obs,name=name,url=img))
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--target", help="one source_name in --targets")
    parser.add_argument("--targets", type=Path)
    parser.add_argument("--deep-pole",action="store_true",help="use two 64-frame NEP windows for a specified target")
    parser.add_argument("--out", type=Path, default=BASE)
    args=parser.parse_args()
    outdir=args.out; outdir.mkdir(parents=True,exist_ok=True)
    ra, dec, pmx, pmy, mag = RA, DEC, PMX, PMY, 16.477
    windows,count,quality,grid,half_width,labels = WINDOWS,64,10,42,38,["Jan 2014","Dec 2023"]
    if args.target:
        target=pd.read_csv(args.targets).set_index("source_name").loc[args.target]
        ra,dec,pmx,pmy,mag=(float(target[k]) for k in ["ra2015_deg","dec2015_deg","pmra","pmdec","w2"])
        windows,count,quality,grid,half_width,labels=[(57065.,57185.),(60352.,60472.)],16,5,72,59,["Spring 2015","Spring 2024"]
    if args.deep_pole:
        if not args.target: raise ValueError("--deep-pole requires --target")
        windows,count,quality,grid,half_width,labels=WINDOWS,64,5,42,38,["Jan 2014","Dec 2023"]
    cd=np.cos(np.deg2rad(dec))
    args.cache.mkdir(parents=True,exist_ok=True)
    manifest=outdir/"manifest.csv"
    rows=list(csv.DictReader(manifest.open())) if manifest.exists() else build_manifest(ra,dec,windows,count,quality)
    east=np.arange(grid,-grid-1,-1)*SCALE
    north=np.arange(-grid+2,grid-1)*SCALE
    X,Y=np.meshgrid(east,north)
    def get(row):
        path=args.cache/(row["name"]+".gz")
        if not path.exists():
            cutout(row["url"],ra,dec,half_width=half_width).writeto(path)
        with fits.open(path) as hdus:
            h=hdus[0].header
            if abs(h["MJD_OBS"]-float(row["mjd"]))>1e-5: raise ValueError("Epoch mismatch")
            for name,key in [("fits_header_sha256","INHDRSHA"),("fits_row_strip_sha256","INRAWSHA")]:
                if row.get(name) and row[name]!=h[key]:raise ValueError("Archive hash changed")
                row[name]=h[key]
            px,py=WCS(h).celestial.world_to_pixel_values(ra+X/(3600*cd),dec+Y/3600)
            out=map_coordinates(np.asarray(hdus[0].data,float),[py,px],order=1,cval=np.nan)
            return (out-np.nanmedian(out))*10**(.4*(20-h["MAGZP"]))
    with ThreadPoolExecutor(8) as pool: images=list(pool.map(get,rows))
    if not manifest.exists():pd.DataFrame(rows).to_csv(manifest,index=False)
    stacks=[]; pos=[]; halves=[]
    for epoch in (0,1):
        indices=[j for j,r in enumerate(rows) if int(r["epoch"])==epoch]
        if len(indices)<4: raise ValueError("Too few frames for both exposure splits")
        g=np.asarray([images[j] for j in indices])
        # Sigma-clipped means retain faint real flux while rejecting particle hits.
        def stack(data):
            med=np.nanmedian(data,axis=0);err=1.4826*np.nanmedian(np.abs(data-med),axis=0)
            return np.nanmean(np.where(np.abs(data-med)<=4*np.maximum(err,1e-8),data,np.nan),axis=0)
        stacks.append(stack(g))
        halves.append([stack(g[i::2]) for i in (0,1)])
        t=np.mean([float(rows[j]["mjd"]) for j in indices])
        pos.append(np.asarray([pmx,pmy])*(t-57170)/365.25)
    # At fixed sky coordinates, a moving source leaves a positive/negative pair in late minus early.
    # Gaussian approximation to the W2 PSF including reprojection; amplitudes in common zero-point DN.
    # Difference local background pedestal and gradient and fit template for the expected dipole.
    sigma=3.25
    def flux(D,ox,oy):
        sources=[]
        support=np.zeros_like(X,dtype=bool)
        for u,v in pos:
            r2=(X-u-ox)**2+(Y-v-oy)**2
            sources.append(np.exp(-r2/(2*sigma**2))*2.75**2/(2*np.pi*sigma**2))
            support|=r2<12.**2
        if np.any(~np.isfinite(D[support])):return np.nan
        A=np.column_stack([(sources[1]-sources[0])[support],np.ones(np.sum(support)),
                           X[support]-ox,Y[support]-oy])
        return np.linalg.lstsq(A,D[support],rcond=None)[0][0]
    D=stacks[1]-stacks[0]
    measured=flux(D,0.,0.)
    controls=[]
    for ox in range(-24,25,8):
        for oy in range(-24,25,8):
            if np.hypot(ox,oy)>=12:controls.append(dict(east=ox,north=oy,flux=float(flux(D,ox,oy))))
    controls=pd.DataFrame(controls);controls.to_csv(outdir/"translations.csv",index=False)
    med=np.nanmedian(controls.flux);mad=1.4826*np.nanmedian(np.abs(controls.flux-med))
    nominal=10**(.4*(20.-mag))
    summary=dict(n_early=sum(int(r["epoch"])==0 for r in rows),n_late=sum(int(r["epoch"])==1 for r in rows),
                 positions_arcsec=[p.tolist() for p in pos],sigma_psf_arcsec=sigma,flux_trace=float(measured),
                 translated_median=float(med),translated_mad=float(mad),nominal_flux=float(nominal),
                 n_controls=len(controls),nominal_minus_observed_over_mad=float((nominal-measured)/mad),
                 split_flux=[float(flux(halves[1][i]-halves[0][i],0.,0.)) for i in (0,1)])
    if not args.target:
        # Use a field source measured by AllWISE to convert approximate Gaussian
        # flux response into an empirical expectation at the candidate magnitude.
        field=pd.read_csv(BASE/"allwise.csv").set_index("designation")
        ref=field.loc["J182533.39+671905.7"]
        rx,ry=(ref.ra-ra)*cd*3600,(ref.dec-dec)*3600
        r2=(X-rx)**2+(Y-ry)**2
        support=r2<12**2
        gauss=np.exp(-r2/(2*sigma**2))*2.75**2/(2*np.pi*sigma**2)
        A=np.column_stack([gauss[support],np.ones(support.sum()),X[support]-rx,Y[support]-ry])
        ref_flux=float(np.linalg.lstsq(A,stacks[1][support],rcond=None)[0][0])
        expected=ref_flux*10**(-.4*(mag-ref.w2mpro))
        summary["field_calibration"]=dict(reference=ref.name,w2=float(ref.w2mpro),
               x_arcsec=rx,y_arcsec=ry,late_reference_flux=ref_flux,
               predicted_dipole_flux=float(expected),gap_over_translation_mad=float((expected-measured)/mad))
    (outdir/"track-test.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    fig,axs=plt.subplots(1,3,figsize=(11,4))
    for i,(ax,im,title) in enumerate(zip(axs,[*stacks,D],[*labels,"Late minus early"])):
        low,high=np.nanpercentile(im,[8,99.4]) if i<2 else (-5,5)
        ax.imshow(im,origin="lower",extent=[east[0],east[-1],north[0],north[-1]],cmap="gray",vmin=low,vmax=high)
        for h,p in enumerate(pos):ax.add_patch(plt.Circle(p,5.5,fill=False,ls=["--","-"][h],color=["#44c2fc","#ff9b30"][h]))
        ax.set_title(title);ax.set_xlabel("East offset (arcsec)")
    axs[0].set_ylabel("North offset (arcsec)")
    fig.tight_layout();fig.savefig(outdir/"difference.png",dpi=130)
    print(json.dumps(summary))


if __name__=="__main__":main()
