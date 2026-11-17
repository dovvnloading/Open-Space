"""Independent image checks of warm literature targets passing the DR1 continuum comparison."""
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"imaging"))
from atlas import BANDS,cloud_cutout,fit
from neighborhoods import cutout, get
sys.path.insert(0,str(ROOT.parent/"wise"))
from fetch import position,separation
from run import read_rows,write_csv,write_json,retry

OUT=ROOT/"scenes"

def wise_scene(task):
    target,band=task
    coadd=target["unwise_objid"].split("o")[0]+"_ac51"
    name=f"{coadd}-w{band}-int-3.fits"
    key=f"{coadd[:2]}/{coadd[:4]}/{coadd}/{name}"
    ra,dec=position(target,2010.5)
    url="https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/allwise/images/p3am_cdd/"+key
    receipt=dict(label=target["label"],band=band,url=url)
    try:
        data,native_hash=retry(cloud_cutout,url,ra,dec)
        receipt["native_sha256"]=native_hash
    except Exception as exc:
        receipt["cloud_error"]=repr(exc)
        url="https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/"+key+"?"+urlencode(
            dict(center=f"{ra:.7f},{dec:.7f}",size="135pix"))
        receipt["fallback_url"]=url
        with urlopen(url,timeout=90) as response: data=response.read()
        if data.startswith(b"\x1f\x8b"): data=gzip.decompress(data)
    with fits.open(io.BytesIO(data)) as hdus:
        arr=hdus[0].data.astype(float)
        w=WCS(hdus[0].header).celestial
    x,y=map(float,w.world_to_pixel_values(ra,dec))
    scale=abs(np.linalg.det(w.pixel_scale_matrix))**.5*3600
    if arr.shape!=(135,135) or min(x,y,134-x,134-y)<16:
        raise ValueError((arr.shape,x,y))
    receipt["cutout_sha256"]=hashlib.sha256(data).hexdigest()
    fname=f"{target['label']}_W{band}.fits.gz"
    (OUT/fname).write_bytes(gzip.compress(data,mtime=0))
    amp,dx,dy,*_=fit(arr,x,y,scale,BANDS[band])
    wide=fit(arr,x,y,scale,BANDS[band]*1.2)
    lon,lat=w.pixel_to_world_values(x+dx,y+dy)
    east=((lon-ra+180)%360-180)*math.cos(math.radians(dec))*3600
    north=(lat-dec)*3600
    yy,xx=np.indices(arr.shape)
    distances=[]
    for ux in (-42,-21,0,21,42):
        for uy in (-42,-21,0,21,42):
            if not(ux or uy): continue
            tx,ty=x+ux,y+uy
            kernel=np.exp(-((xx-tx)**2+(yy-ty)**2)/(2*(BANDS[band]/scale/2.355)**2))
            offsets=fit(arr+amp*kernel,tx,ty,scale,BANDS[band])[1:3]
            distances.append(scale*np.hypot(*offsets))
    row=dict(label=target["label"],band=band,east_arcsec=round(float(east),3),north_arcsec=round(float(north),3),
             gaia_offset_arcsec=round(float(np.hypot(east,north)),3),
             wider_psf_change_arcsec=round(float(scale*np.hypot(wide[1]-dx,wide[2]-dy)),3),
             blank_placements=len(distances),placement_p68_arcsec=round(float(np.percentile(distances,68)),3),
             placement_p90_arcsec=round(float(np.percentile(distances,90)),3),
             at_fit_bound=max(abs(float(dx)),abs(float(dy)))*scale>=6.98)
    return row,receipt,(arr,x,y,scale,dx,dy)


def display_scene(target,layer):
    rec,img=cutout((target,layer))
    if layer=="ls-dr10" and rec["status"]!="displayed":
        url="https://alasky.cds.unistra.fr/hips-image-services/hips2fits?"+urlencode(dict(
            hips="CDS/P/DSS2/red",width=120,height=120,ra=target["ra"],dec=target["dec"],
            fov=48/3600,projection="TAN",format="fits"))
        try:
            with urlopen(url,timeout=50) as response: raw=response.read()
            with fits.open(io.BytesIO(raw)) as hdus: a=hdus[0].data.astype(float).squeeze()
            lo,hi=np.nanpercentile(a,[3,99])
            if a.shape!=(120,120) or not hi>lo: raise ValueError("empty DSS image")
            img=Image.fromarray(np.uint8(np.clip((a-lo)/(hi-lo),0,1)*255)[::-1]).convert("RGB")
            rec.update(dss_fallback_url=url,dss_sha256=hashlib.sha256(raw).hexdigest(),optical_display="DSS2 red")
        except Exception as exc: rec["dss_error"]=repr(exc)
    return rec,img


def gaia_neighbors(targets):
    cones=[f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r['ra']},{r['dec']},{20/3600:.9f}))" for r in targets]
    query="SELECT source_id,ra,dec,phot_g_mean_mag,parallax,pmra,pmdec,ruwe,non_single_star FROM gaiadr3.gaia_source WHERE "+" OR ".join(cones)
    url="https://gea.esac.esa.int/tap-server/tap/sync?"+urlencode(dict(LANG="ADQL",FORMAT="csv",REQUEST="doQuery",QUERY=query))
    with urlopen(url,timeout=90) as resp: data=resp.read()
    catalog=list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
    candidates=[]
    for r in targets:
        for s in catalog:
            arc=separation(float(s["ra"]),float(s["dec"]),float(r["ra"]),float(r["dec"]))
            if arc<=20.01:
                candidates.append(dict(label=r["label"],offset_arcsec=round(arc,3),**s))
    if not all(any(c["source_id"]==r["source_id"] for c in candidates) for r in targets): raise ValueError("Gaia targets missing")
    write_csv(OUT/"gaia_neighbors.csv",sorted(candidates,key=lambda c:(c["label"],c["offset_arcsec"])))
    return dict(url=url,sha256=hashlib.sha256(data).hexdigest(),rows=len(catalog))


def main():
    OUT.mkdir(exist_ok=True)
    astrom={r["source_id"]:r for r in read_rows(ROOT/"positions.csv")}
    baseline={r["source_id"]:r for r in read_rows(ROOT.parent/"literature/fgk_baseline.csv")}
    stars=[dict(**astrom[r["source_id"]],unwise_objid=baseline[r["source_id"]]["unwise_objid"])
           for r in read_rows(ROOT/"all_53.csv") if 100<=float(r["published_bb_k"])<=600
           and r["continuum_check"]=="redder_than_both_in_both_masks"]
    if not stars: raise ValueError("No warm flagged fields")
    if "--plot-only" in sys.argv:
        fit_rows = {(r["label"], int(r["band"])): r for r in read_rows(OUT/"centroids.csv")}
        wise = []
        for target in stars:
            ra, dec = position(target, 2010.5)
            for band in BANDS:
                label = target["label"]
                with fits.open(OUT/f"{label}_W{band}.fits.gz") as hdus:
                    arr = hdus[0].data.astype(float)
                    w = WCS(hdus[0].header).celestial
                x, y = map(float, w.world_to_pixel_values(ra, dec))
                scale = abs(np.linalg.det(w.pixel_scale_matrix))**.5 * 3600
                row = fit_rows[label, band]
                east = float(row["east_arcsec"]) / (3600 * math.cos(math.radians(dec)))
                north = float(row["north_arcsec"]) / 3600
                cx, cy = map(float, w.world_to_pixel_values(ra + east, dec + north))
                wise.append((None, None, (arr, x, y, scale, cx-x, cy-y)))
        display = []
        for rec in json.loads((OUT/"viewer_queries.json").read_text()):
            if rec["status"] == "displayed":
                raw = get(rec["url"])
                if hashlib.sha256(raw).hexdigest() != rec["sha256"]:
                    raise ValueError(f"Viewer image changed for {rec['label']} {rec['layer']}")
                img = Image.open(io.BytesIO(raw)).convert("RGB")
            elif "dss_sha256" in rec:
                raw = get(rec["dss_fallback_url"])
                if hashlib.sha256(raw).hexdigest() != rec["dss_sha256"]:
                    raise ValueError(f"DSS image changed for {rec['label']}")
                with fits.open(io.BytesIO(raw)) as hdus:
                    arr = hdus[0].data.astype(float).squeeze()
                lo, hi = np.nanpercentile(arr, [3, 99])
                img = Image.fromarray(np.uint8(np.clip((arr-lo)/(hi-lo), 0, 1)*255)[::-1]).convert("RGB")
            else:
                img = Image.new("RGB", (120, 120), "#202020")
            display.append((rec, img))
    else:
        grouped=[(t,b) for t in stars for b in BANDS]
        with ThreadPoolExecutor(4) as pool: wise=list(pool.map(wise_scene,grouped))
        write_csv(OUT/"centroids.csv",[r[0] for r in wise])
        write_json(OUT/"wise_queries.json",[r[1] for r in wise])
        write_json(OUT/"gaia_query.json",retry(gaia_neighbors,stars))
        cases=[(t,l) for t in stars for l in ("ls-dr10","unwise-neo6")]
        with ThreadPoolExecutor(4) as pool: display=list(pool.map(lambda t:display_scene(*t),cases))
        write_json(OUT/"viewer_queries.json",[r[0] for r in display])
    fig,axes=plt.subplots(len(stars),5,figsize=(12.6,2.65*len(stars)+.8),squeeze=False)
    for idx,t in enumerate(stars):
        for col,layer in enumerate(("Legacy Surveys DR10","unWISE")):
            rec,img=display[2*idx+col]
            ax=axes[idx,col]; ax.imshow(img); ax.set(xticks=[],yticks=[])
            title=rec.get("optical_display",layer)
            if idx == 0: ax.set_title(layer, fontsize=9)
            if rec.get("status") != "displayed" and "optical_display" not in rec:
                ax.text(.5, .42, "Unavailable", color="white", transform=ax.transAxes,
                        ha="center", va="center", fontsize=9)
            elif title != layer:
                ax.text(.03, .04, title, color="white", transform=ax.transAxes, fontsize=9)
            ax.plot(59.5,59.5,marker="+",color="cyan",ms=9)
            if col==0: ax.set_ylabel(t["label"])
        for j,band in enumerate(BANDS):
            arr,x,y,s,dx,dy=wise[3*idx+j][2]; ax=axes[idx,2+j]
            p1,p2=np.nanpercentile(arr,[15,99.2])
            ax.imshow(np.arcsinh((arr-p1)*4/max(p2-p1,1e-6)),origin="lower",cmap="gray",vmin=-.1,vmax=2.6)
            ax.plot(x,y,marker="+",ms=10,color="cyan");ax.plot(x+dx,y+dy,marker="o",fillstyle="none",ms=8,color="orange")
            ax.set(xlim=(x-24/s,x+24/s),ylim=(y-24/s,y+24/s),xticks=[],yticks=[])
            if idx == 0: ax.set_title(f"AllWISE W{band}", fontsize=9)
    fig.suptitle("Published warm targets with a red DR1 continuum",fontsize=15)
    fig.text(.5,.015,"48 arcsec fields. Cyan marks Gaia. Orange marks the fitted AllWISE centroid.",ha="center",fontsize=9)
    fig.subplots_adjust(left=.038, right=.992, bottom=.067, top=.865, wspace=.055, hspace=.07)
    fig.savefig(OUT/"fields.png",dpi=130);plt.close(fig)

if __name__=="__main__": main()
