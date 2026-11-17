"""Forced mass-centre photometry on calibrated Q1 cutouts and published MGE source models."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from cutouts import CutoutArchive, image_extensions

BANDS = ("VIS","NIR_Y","NIR_J","NIR_H")
HERE = Path(__file__).resolve().parent


def aperture(radius, y, x):
    """Same whole-pixel circular aperture at object and blank positions."""
    r = int(np.ceil(radius))+1
    yy,xx = np.mgrid[-r:r+1,-r:r+1]
    return (yy[(yy-(y-round(y)))**2+(xx-(x-round(x)))**2 <= radius**2],
            xx[(yy-(y-round(y)))**2+(xx-(x-round(x)))**2 <= radius**2])


def blank_measurement(image, rms, valid, y, x, radius):
    dy, dx = aperture(radius, y, x)
    step = max(9,int(np.ceil(2*radius))+2)
    sky = []
    for cy in range(step, image.shape[0]-step, step):
        for cx in range(step, image.shape[1]-step, step):
            dist = np.hypot(cy-y,cx-x)
            if dist < 70+radius or dist > 135-radius:
                continue
            yy,xx=cy+dy,cx+dx
            if valid[yy,xx].all():
                sky.append(image[yy,xx].sum(dtype=float))
    if len(sky)<12:
        return dy,dx,np.nan,np.nan,len(sky)
    median=np.median(sky)
    sigma=1.4826022*np.median(np.abs(np.asarray(sky)-median))
    return dy,dx,median,sigma,len(sky)


def ab(zero, f):
    return zero-2.5*np.log10(f) if f>0 else np.nan


def stamp(draws, obj, out):
    """Rendering only: each stamp has a separate asinh display stretch."""
    board=Image.new("RGB",(176*len(draws),196),"white")
    d=ImageDraw.Draw(board)
    for i,(name,arr) in enumerate(draws):
        if arr is None:
            d.text((176*i+8,88),"missing model",fill="black")
        else:
            lo,hi=np.nanpercentile(arr,[10,99.3])
            stretch=np.arcsinh(10*np.maximum(arr-lo,0)/max(hi-lo,1e-8))/np.arcsinh(10)
            stampim=Image.fromarray(np.uint8(255*np.nan_to_num(stretch).clip(0,1))[::-1]).resize((176,176))
            board.paste(stampim,(176*i,0))
        d.text((176*i+4,177),name,fill="black")
    board.save(out/f"{obj}.png",optimize=True)


def analyze(obj, mass, archive, make_stamp, out):
    names = [obj+".fits","info.json","result/mge_lens_light.fits","result/mge_source_light.fits"]
    data = archive.members(obj,names)
    hashes = [{"id_str":obj,"archive":archive.name+".zip","member":f"{archive.name}/{obj}/{k}",
               "bytes":len(v),"sha256":hashlib.sha256(v).hexdigest()} for k,v in data.items()]
    raw = image_extensions(data[obj+".fits"])
    lens = image_extensions(data[names[2]]) if names[2] in data else {}
    src = image_extensions(data[names[3]]) if names[3] in data else {}
    info = json.loads(data["info.json"])
    centre_y=float(mass["centre_0_median_pdf"])
    centre_x=float(mass["centre_1_median_pdf"])
    theta=float(mass["einstein_radius_effective_median_pdf"])
    if not np.isfinite([theta,centre_y,centre_x]).all():
        raise ValueError("non-finite mass centre or theta")
    rows,draws=[],{}
    for band in BANDS:
        if band+"_FLUX" not in raw:
            continue
        h,im=raw[band+"_FLUX"]
        im=im.astype(float)
        rms=raw[band+"_RMS"][1]
        valid=np.isfinite(im)&np.isfinite(rms)&(rms>0)&(rms<1e15)
        pix=float(abs(float(h["CD2_2"]))*3600)
        ny,nx=im.shape
        # The release pipeline crops the central 150x150 pixels. Its FITS writer
        # uses increasing y with increasing FITS row, increasing x with column.
        y,x=(ny-1)/2+centre_y/pix,(nx-1)/2+centre_x/pix
        if not (0<=y<ny and 0<=x<nx):
            raise ValueError("Published centre outside cutout")
        model=src.get(band,(None,None))[1]
        light=lens.get(band,(None,None))[1]
        def aligned(m):
            if m is None: return None
            oy,ox=(ny-m.shape[0])//2,(nx-m.shape[1])//2
            if (ny-m.shape[0])%2 or (nx-m.shape[1])%2 or min(oy,ox)<0:
                raise ValueError("Unrecognized model dimensions")
            canvas=np.full(im.shape,np.nan)
            canvas[oy:oy+m.shape[0],ox:ox+m.shape[1]]=m
            return canvas
        model,light=aligned(model),aligned(light)
        window=(slice(round(y)-35,round(y)+36),slice(round(x)-35,round(x)+36))
        draws[(band,"raw")]=im[window]
        draws[(band,"residual")]=None if light is None else (im-light)[window]
        for label,a in [("core",0.3),("einstein",theta)]:
            if a<=0 or a/pix>min(x,nx-1-x,y,ny-1-y): continue
            dy,dx,bg,sky_sigma,nsky=blank_measurement(im,rms,valid,y,x,a/pix)
            yy,xx=round(y)+dy,round(x)+dx
            counts=im[yy,xx].sum()-bg if valid[yy,xx].all() else np.nan
            src_flux=np.sum(model[yy,xx]) if model is not None else np.nan
            lens_flux=np.sum(light[yy,xx]) if light is not None else np.nan
            err_diag=np.sqrt(np.sum(np.square(rms[yy,xx].astype(float))))
            good=bool(valid[yy,xx].all())
            sigma=max(sky_sigma,err_diag) if good and np.isfinite(sky_sigma) else np.nan
            zero=float(h["MAGZERO"])
            rows.append(dict(id_str=obj,ra_cutout=float(info["target_ra"]),dec_cutout=float(info["target_dec"]),
                band=band,aperture=label,radius_arcsec=a,mass_dy_arcsec=centre_y,mass_dx_arcsec=centre_x,
                npix=len(dx),coverage_ok=good,n_sky=nsky,zp_ab=zero,aperture_sky_counts=bg,aperture_sky_mad=sky_sigma,rms_quadrature=err_diag,
                observed_counts=counts,source_model_counts=src_flux,lens_model_counts=lens_flux,
                source_subtracted_counts=counts-src_flux,observed_snr=counts/sigma,
                source_subtracted_snr=(counts-src_flux)/sigma,source_subtracted_mag=ab(zero,counts-src_flux),
                lens_model_mag=ab(zero,lens_flux),noise_only_5sigma_mag=ab(zero,5*sigma)))
    if make_stamp:
        stamp([(b+" observed" if k=="raw" else b+" minus lens",draws.get((b,k)))
               for b,k in [("VIS","raw"),("VIS","residual"),("NIR_H","raw"),("NIR_H","residual")]],obj,out)
    return rows,hashes


def table(path, rows):
    rows.sort(key=lambda x:tuple(str(x.get(k,"")) for k in ("id_str","band","aperture","member")))
    temp=path.with_suffix(".tmp")
    with temp.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader()
        for r in rows:
            w.writerow({k:format(v,".9g") if isinstance(v,(float,np.floating)) and np.isfinite(v) else
                        ("" if isinstance(v,(float,np.floating)) else v) for k,v in r.items()})
    temp.replace(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mass_csv",type=Path)
    p.add_argument("--index-cache",type=Path)
    p.add_argument("--out",type=Path,default=HERE)
    p.add_argument("--all",action="store_true")
    p.add_argument("--workers",type=int,default=4)
    args=p.parse_args()
    out=args.out
    out.mkdir(parents=True,exist_ok=True)
    priority=set((HERE/"priority.txt").read_text().split())
    mass={r["id_str"]:r for r in csv.DictReader(args.mass_csv.open())}
    stamps=priority | {"102158887_2655306065654633189","102042915_NEG528545016288651034"}
    ids=[i for i in mass if (i in priority or args.all) and mass[i]["centre_0_median_pdf"] and mass[i]["einstein_radius_effective_median_pdf"]]
    if not priority.issubset(ids):
        raise ValueError("Missing a priority ID")
    arc=CutoutArchive("lens",args.index_cache)
    rows,hashes=[],[]
    if args.all and (out/"apertures.csv").exists() and "coverage_ok" in (out/"apertures.csv").open().readline() and (out/"input_members.csv").exists():
        rows=list(csv.DictReader((out/"apertures.csv").open()))
        hashes=list(csv.DictReader((out/"input_members.csv").open()))
        done={r["id_str"] for r in rows} & {r["id_str"] for r in hashes}
        ids=[i for i in ids if i not in done]
        rows=[r for r in rows if r["id_str"] in done]
        hashes=[r for r in hashes if r["id_str"] in done]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(analyze,i,mass[i],arc,i in stamps,out):i for i in ids}
        errors=[]
        for j,f in enumerate(as_completed(futures),1):
            try:
                r,h=f.result(); rows.extend(r);hashes.extend(h)
            except Exception as e:
                errors.append({"id_str":futures[f],"error":repr(e)})
            if (j%10==0 or j==len(ids)) and rows:
                table(out/"apertures.csv",rows);table(out/"input_members.csv",hashes)
                print(f"{j}/{len(ids)} attempted, {len(errors)} errors",flush=True)
    if errors:
        table(out/"failed.csv",errors)
        raise SystemExit(f"{len(errors)} objects failed")


if __name__=="__main__":
    main()
