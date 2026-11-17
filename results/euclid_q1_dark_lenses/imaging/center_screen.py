"""Measure blended light at the released cutout origin for A/B systems without a finite mass-centre fit."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from measure import BANDS,HERE,ab,blank_measurement,table
from cutouts import CutoutArchive,image_extensions


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("catalog",type=Path)
    p.add_argument("mass",type=Path)
    p.add_argument("--index-cache",type=Path)
    args=p.parse_args()
    targets={r["id_str"] for r in csv.DictReader(args.catalog.open()) if r["subset"]=="discovery_engine" and r["grade"] in ("A","B")}
    fits={r["id_str"] for r in csv.DictReader(args.mass.open()) if r["einstein_radius_effective_median_pdf"] and r["centre_0_median_pdf"]}
    need=targets-fits
    catalogues={}
    for name in ("lens","recenter","group","unsuccess"):
        arc=CutoutArchive(name,args.index_cache)
        for id_ in need:
            if f"{name}/{id_}/{id_}.fits" in arc.entries and id_ not in catalogues:catalogues[id_]=arc
    def check(item):
        obj,arc=item
        name=obj+".fits"
        d=arc.members(obj,[name,"info.json"])
        raw=image_extensions(d[name])
        try:
            info=json.loads(d.get("info.json",b"{}"))
        except (ValueError, UnicodeError):
            info={}
        rows=[]
        for b in BANDS:
            if b+"_FLUX" not in raw:continue
            h,im=raw[b+"_FLUX"]
            rms=raw[b+"_RMS"][1]
            y,x=(np.array(im.shape)-1)/2
            dy,dx,bg,scat,n=blank_measurement(im,rms,np.isfinite(im)&np.isfinite(rms)&(rms>0)&(rms<1e15),y,x,0.3/(3600*abs(float(h["CD2_2"]))))
            yy,xx=round(y)+dy,round(x)+dx
            f=np.sum(im[yy,xx],dtype=float)-bg
            if not np.all(np.isfinite(rms[yy,xx])&(rms[yy,xx]>0)&(rms[yy,xx]<1e15)):f=np.nan
            sigma=max(scat,np.sqrt(np.sum(rms[yy,xx].astype(float)**2)))
            rows.append(dict(id_str=obj,archive=arc.name,ra_cutout=info.get("target_ra"),dec_cutout=info.get("target_dec"),
                             band=b,sky_apertures=n,blended_counts=f,sigma_counts=sigma,blended_snr=f/sigma,
                             blended_mag=ab(float(h["MAGZERO"]),f),noise_only_5sigma_mag=ab(float(h["MAGZERO"]),5*sigma) if np.isfinite(f) else np.nan))
        hh=[dict(id_str=obj,archive=arc.name+".zip",member=f"{arc.name}/{obj}/{k}",bytes=len(v),sha256=hashlib.sha256(v).hexdigest()) for k,v in d.items()]
        return rows,hh
    rows,hashes=[],[]
    failures=[dict(id_str=i,reason="archive cutout unavailable") for i in sorted(need-catalogues.keys())]
    with ThreadPoolExecutor(5) as pool:
        tasks={pool.submit(check,item):item[0] for item in sorted(catalogues.items())}
        for j,task in enumerate(as_completed(tasks),1):
            try:
                r,h=task.result()
                rows.extend(r);hashes.extend(h)
                if not r:failures.append(dict(id_str=tasks[task],reason="zero-byte FITS or no calibrated flux planes"))
            except (ValueError,KeyError,OSError) as e:
                failures.append(dict(id_str=tasks[task],reason=str(e)))
            if j%20==0 or j==len(catalogues):
                table(HERE/"unmodelled_centres.csv",rows);table(HERE/"unmodelled_inputs.csv",hashes)
                print(j,"/",len(catalogues),flush=True)
    if failures:table(HERE/"unmodelled_unavailable.csv",failures)
    print("unavailable",len(failures))


if __name__=="__main__":main()
