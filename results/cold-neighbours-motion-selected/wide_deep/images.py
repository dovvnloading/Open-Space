"""Test a surviving CatWISE motion in early and held-out NEOWISE single exposures."""
import argparse,csv,hashlib,math,os,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.ndimage import map_coordinates
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent/"pilot"/"j0352-images"))
from native import cutout

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source");parser.add_argument("--catalog",type=Path,default=BASE.parent/"wide"/"selected.csv")
    parser.add_argument("--cache",type=Path,default=BASE/"native-cache")
    a=parser.parse_args();s=next(r for r in csv.DictReader(a.catalog.open()) if r["source_name"]==a.source)
    ra,dec,ma,md=map(float,(s["ra2015_deg"],s["dec2015_deg"],s["pmra"],s["pmdec"]));mag=float(s["w2"])
    out=BASE/a.source;out.mkdir(exist_ok=True);cache=a.cache/a.source;cache.mkdir(parents=True,exist_ok=True)
    url="https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?"+urlencode(dict(POS=f"{ra},{dec}",SIZE=".001",COLUMNS="scan_id,frame_num,band,qual_frame,mjd_obs,crval1,crval2"))
    (out/"frames.url").write_text(url+"\n");path=out/"frames.tbl"
    if not path.exists():path.write_bytes(urlopen(url,timeout=120).read())
    t=Table.read(path,format="ascii.ipac");t=t[(t["band"]==2)&(t["qual_frame"]>=5)];t.sort("mjd_obs")
    visits=np.split(np.arange(len(t)),np.flatnonzero(np.diff(t["mjd_obs"])>50)+1)
    if len(visits)<12:raise ValueError("Fewer than 12 visits")
    med=[np.median(t["mjd_obs"][v]) for v in visits];cos=math.cos(math.radians(dec));rows=[]
    ev=[k for k,m in enumerate(med) if m<57400][:4]
    lv=[k for k,m in enumerate(med) if m>59215][-4:]
    for k in ev+lv:
        block=t[visits[k]];sep=np.hypot((block["crval1"]-ra)*cos,block["crval2"]-dec)
        for f in block[np.lexsort((sep,-np.asarray(block["qual_frame"])))[:8]]:
            scan=f["scan_id"];frame=f'{f["frame_num"]:03d}';fn=f"{scan}{frame}-w2-int-1b.fits"
            rows.append(dict(visit=k,mjd=float(f["mjd_obs"]),qual_frame=int(f["qual_frame"]),file=fn,url=f"https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/neowiser/images/p1bm_frm/{scan[-2:]}/{scan}/{frame}/{fn}"))
    ax=np.arange(-40,40.01,1.375);xx,yy=np.meshgrid(ax,ax);fw=6.4
    def get(r):
        p=cache/(r["file"]+".gz")
        if not p.exists():cutout(r["url"],ra,dec).writeto(p)
        with fits.open(p) as h:h0=h[0].header;im=h[0].data.astype(float)
        if abs(h0["MJD_OBS"]-r["mjd"])>1e-5:raise ValueError("wrong timestamp")
        r.update(header_sha256=h0["INHDRSHA"],row_strip_sha256=h0["INRAWSHA"])
        w=WCS(h0).celestial;u,v=w.world_to_pixel_values(ra+xx/(cos*3600),dec+yy/3600);dt=(r["mjd"]-57170)/365.25
        cx,cy=w.world_to_pixel_values(ra+ma*dt/(cos*3600),dec+md*dt/3600)
        gy,gx=np.indices(im.shape);sx,sy=proj_plane_pixel_scales(w)*3600
        model=np.exp(-4*np.log(2)*(((gx-cx)*sx)**2+((gy-cy)*sy)**2)/fw**2)
        both=[]
        for ar in [im,im+model/model.sum()*10**(.4*(h0["MAGZP"]-mag))]:
            proj=map_coordinates(ar,[v,u],order=1,cval=np.nan);both.append((proj-np.nanmedian(proj))*10**(.4*(20.-h0["MAGZP"])))
        return r,both
    with ThreadPoolExecutor(6) as ex:dat=list(ex.map(get,rows))
    def write(name,rs):
        with (out/name).open("w") as f:
            w=csv.DictWriter(f,fieldnames=rs[0]);w.writeheader();w.writerows(rs)
    write("inputs.csv",rows);fig,axs=plt.subplots(2,3,figsize=(8,5.5),sharex=True,sharey=True);scores=[]
    for season in range(2):
        eg=[k for k in ev if k%2==season];lg=[k for k in lv if k%2==season]
        if not eg or not lg:raise ValueError("scan seasons missing")
        dates=[np.median([r["mjd"] for r in rows if r["visit"] in z]) for z in (eg,lg)]
        dt1,dt2=(np.array(dates)-57170)/365.25;xe,ye=ma*dt1,md*dt1;xl,yl=ma*dt2,md*dt2
        def contrast(d,ox,oy):
            g=lambda x,y:np.exp(-4*np.log(2)*((xx-x-ox)**2+(yy-y-oy)**2)/fw**2)
            model=g(xl,yl)-g(xe,ye);return np.nansum(d*model)/np.nansum(model**2)
        vals=[]
        for ix in [0,1]:
            A,B=[np.nanmedian([v[ix] for r,v in dat if r["visit"] in block],axis=0) for block in (eg,lg)]
            D=B-A;vals.append(contrast(D,0,0))
            if ix==1:continue
            ctrl=[contrast(D,ox,oy) for ox in range(-25,26,7) for oy in range(-25,26,7) if np.hypot(ox,oy)>=12]
            medc=np.median(ctrl);mad=1.4826*np.median(np.abs(ctrl-medc))
            for j,im in enumerate((A,B,D)):
                off=im[np.hypot(xx,yy)>10];me=np.nanmedian(off);si=1.4826*np.nanmedian(np.abs(off-me))
                panel=axs[season,j];panel.imshow(im,extent=[ax.min()-.6875,ax.max()+.6875,ax.min()-.6875,ax.max()+.6875],origin="lower",cmap="gray",vmin=me-2*si,vmax=me+5*si)
                for x,y,c in [(xe,ye,"#eab308"),(xl,yl,"#5eead4")]:panel.add_patch(plt.Circle((x,y),3.2,fill=False,lw=1,ec=c))
                panel.set_xlim(26,-26);panel.set_ylim(-26,26)
                if season==0:panel.set_title(["Early visits","Late visits","Late minus early"][j])
                if j==0:panel.set_ylabel(f"Season {season+1}\nNorth (arcsec)")
                if season==1:panel.set_xlabel("East (arcsec)")
        scores.append(dict(season=season,early_mjd=dates[0],late_mjd=dates[1],amplitude=vals[0],off_median=medc,off_mad=mad,raw_over_mad=(vals[0]-medc)/mad,injected_over_mad=(vals[1]-medc)/mad))
    write("test.csv",scores);fig.suptitle(a.source+" · NEOWISE W2");fig.tight_layout();fig.savefig(out/"epochs.png",dpi=140)
    print(a.source,[(round(i["raw_over_mad"],2),round(i["injected_over_mad"],2)) for i in scores])
if __name__=="__main__":main()
