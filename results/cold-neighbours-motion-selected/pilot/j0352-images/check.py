"""Compare scan-matched native W2 stacks at the faint partial-track test location."""
import csv,hashlib,math,os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from scipy.ndimage import map_coordinates
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from native import cutout

HERE=Path(__file__).resolve().parent
RA,DEC,MA,MD=58.0379644,-49.9902745,.030636,-.570816
COS=np.cos(np.deg2rad(DEC))
FWHM=6.4
AX=np.arange(-40,40.01,1.375)
X,Y=np.meshgrid(AX,AX)
CACHE=Path(os.environ.get("J0352_CACHE",HERE/"cache"))
def manifest():
    with (HERE/"frames.url").open() as f:url=f.read().strip()
    if not (HERE/"frames.tbl").exists(): (HERE/"frames.tbl").write_bytes(urlopen(url).read())
    t=Table.read(HERE/"frames.tbl",format="ascii.ipac")
    t=t[(t["band"]==2)&(t["qual_frame"]==10)];t.sort("mjd_obs")
    visits=np.split(np.arange(len(t)),np.flatnonzero(np.diff(t["mjd_obs"])>45)+1)
    rows=[]
    for i in [0,1,2,3,18,19,20,21]:
        g=t[visits[i]]
        sep=np.hypot((g["crval1"]-RA)*COS,g["crval2"]-DEC)
        for r in g[np.argsort(sep)[:16]]:
            scan=r["scan_id"];frame=f'{r["frame_num"]:03d}'
            name=f"{scan}{frame}-w2-int-1b.fits"
            url=("https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/"
                 f"p1bm_frm/{scan[-2:]}/{scan}/{frame}/{name}")
            rows.append(dict(visit=i,mjd=float(r["mjd_obs"]),file=name,url=url))
    return rows
def get(r):
    path=CACHE/(r["file"]+".gz")
    if not path.exists():cutout(r["url"],RA,DEC).writeto(path)
    with fits.open(path) as hdus:h=hdus[0].header;data=hdus[0].data.copy()
    if abs(h["MJD_OBS"]-r["mjd"])>1e-5:raise ValueError("epoch mismatch")
    r.update(header_sha256=h["INHDRSHA"],row_strip_sha256=h["INRAWSHA"])
    u,v=WCS(h).celestial.world_to_pixel_values(RA+X/(COS*3600),DEC+Y/3600)
    arr=map_coordinates(data.astype(float),[v,u],order=1,cval=np.nan)
    return r, (arr-np.nanmedian(arr))*10**(.4*(20-h["MAGZP"]))
def psf(x,y):return np.exp(-4*np.log(2)*((X-x)**2+(Y-y)**2)/FWHM**2)
def amp(d,xa,ya,xb,yb):
    v=psf(xb,yb)-psf(xa,ya)
    return np.nansum(d*v)/np.nansum(v*v)
def main():
    CACHE.mkdir(exist_ok=True,parents=True)
    rows=manifest()
    with ThreadPoolExecutor(8) as p: vals=list(p.map(get,rows))
    with (HERE/"inputs.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    fig,axes=plt.subplots(2,3,figsize=(8.5,5.5),sharex=True,sharey=True)
    results=[]
    for i,(early,late) in enumerate([([0,2],[18,20]),([1,3],[19,21])]):
        imgs=[];dates=[]
        for epochs in [early,late]:
            sub=[(r,a) for r,a in vals if r["visit"] in epochs]
            imgs.append(np.nanmedian([a for r,a in sub],axis=0))
            dates.append(np.median([r["mjd"] for r,a in sub]))
        a,b=imgs;d=b-a;t1,t2=(np.array(dates)-57170)/365.25
        xe,ye=t1*MA,t1*MD;xl,yl=t2*MA,t2*MD
        c=amp(d,xe,ye,xl,yl)
        control=[amp(d,xe+ox,ye+oy,xl+ox,yl+oy) for ox in np.arange(-25,26,7) for oy in np.arange(-25,26,7) if np.hypot(ox,oy)>=12]
        med=np.median(control);mad=1.4826*np.median(np.abs(control-med))
        results.append(dict(group=i,early_mjd=dates[0],late_mjd=dates[1],early_frames=sum(r["visit"] in early for r in rows),late_frames=sum(r["visit"] in late for r in rows),amplitude=c,off_median=med,off_mad=mad,contrast_over_mad=(c-med)/mad,off_positions=len(control),n_off_at_least_target=sum(x>=c for x in control)))
        for j,im in enumerate([a,b,d]):
            ax=axes[i,j];outer=im[np.hypot(X,Y)>10];m=np.nanmedian(outer);s=1.4826*np.nanmedian(np.abs(outer-m))
            ax.imshow(im,extent=[AX.min()-.6875,AX.max()+.6875,AX.min()-.6875,AX.max()+.6875],origin="lower",cmap="gray",vmin=m-2*s,vmax=m+5*s)
            for x,y,col in [(xe,ye,"#eab308"),(xl,yl,"#5eead4")]:
                ax.add_patch(plt.Circle((x,y),3.2,fill=False,lw=1,ec=col))
            if i==0:ax.set_title(["2014–15","2023–24","Late minus early"][j])
            ax.set_xlim(26,-26);ax.set_ylim(-26,26)
            if j==0:ax.set_ylabel(f"Season {i+1}\nNorth (arcsec)")
            if i==1:ax.set_xlabel("East (arcsec)")
    fig.suptitle("J035209.10−495925.4 · NEOWISE W2")
    fig.tight_layout();fig.savefig(HERE/"epochs.png",dpi=140);plt.close(fig)
    with (HERE/"test.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
    print(results)
if __name__=="__main__":main()
