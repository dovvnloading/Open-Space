#!/usr/bin/env python3
"""Render local and full-stamp Euclid VIS views of positive host-matched segments."""
import argparse
import csv
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from PIL import Image,ImageDraw,ImageFont


def render(image, valid, rect, side, points, flip=False):
    x,y,w,h=rect
    raw=image[y:y+h,x:x+w]
    ok=valid[y:y+h,x:x+w]
    a=raw[ok]
    median=np.median(a) if a.size else 0
    noise=np.median(np.abs(a-median))*1.4826 if a.size else 1
    noise=max(noise,1e-9)
    intensity=np.arcsinh(np.where(ok,np.maximum(0,(raw-median)/noise+1),0))
    grey=np.asarray(255*np.clip(intensity/np.arcsinh(25),0,1),dtype=np.uint8)
    if flip: grey=np.flipud(grey)
    rgb=Image.fromarray(grey).resize((side,side),Image.Resampling.LANCZOS).convert("RGB")
    d=ImageDraw.Draw(rgb)
    for (px,py),color,radius in zip(points,["#64edf2","#64edf2","#ffb23c"],[6,6,9]):
        xx=(px-x+.5)/w*side-.5
        yy=(py-y+.5)/h*side-.5
        if flip: yy=side-1-yy
        d.ellipse((xx-radius,yy-radius,xx+radius,yy+radius),outline=color,width=1)
    return rgb


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts",required=True,type=Path)
    p.add_argument("--photometry",required=True,type=Path)
    p.add_argument("--images",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--all-hosts",action="store_true",help="include matches failing the VIS flux gate")
    p.add_argument("--skip-reviewed",type=Path,help="optional CSV with previously inspected IDs")
    a=p.parse_args()
    done=set()
    if a.skip_reviewed:
        with a.skip_reviewed.open() as f:
            done={r["id"] for r in csv.DictReader(f)}
    with a.hosts.open() as f: hosts=list(csv.DictReader(f))
    with a.photometry.open() as f:
        photo={r["id"]:r for r in csv.DictReader(f) if r["band"]=="VIS"}
    targets=[]
    for row in hosts:
        if row["id"] in done: continue
        if row["sign"]!="positive" or not (row["host_object_id"] or row.get("axis_host_id")): continue
        ph=photo.get(row["id"],{})
        if not a.all_hosts and not (ph.get("usable")=="1" and int(ph["positive_bins"])>=4
                                    and float(ph["central_ridge_contrast_sigma"])>=3):
            continue
        targets.append(row)
    targets.sort(key=lambda r:(r["tile"],int(r["index"])))
    a.out.mkdir(exist_ok=True,parents=True)
    side,cols,rows_per_page=250,2,5
    cw,ch=2*side+18,side+49
    paths=[]
    for i,row in enumerate(targets):
        index=i%(cols*rows_per_page)
        if index==0:
            page=Image.new("RGB",(cols*cw,26+ch*rows_per_page),"white")
            draw=ImageDraw.Draw(page)
            draw.text((8,6),"VIS  |  left: full stamp   right: local detail   cyan: endpoints   amber: catalogued host   north up, east left",fill="black")
        xx=(index%cols)*cw
        yy=26+(index//cols)*ch
        tile=row["tile"]
        with fits.open(a.images/f"{tile}_VIS.fits") as vis, fits.open(a.images/f"{tile}_rms.fits") as rms:
            im=vis[0].data
            r=rms[0].data
            valid=np.isfinite(im)&(im!=0)&np.isfinite(r)&(r>0)&(r<1e8)
            wcs=WCS(vis[0].header)
            if row["host_object_id"]:
                host=wcs.world_to_pixel_values(float(row["host_ra"]),float(row["host_dec"]))
                host_label=f'endpoint offset {float(row["host_sep_arcsec"]):.2f}"'
            else:
                # Reconstruct projected host position from the host audit's along/cross
                # coordinates when its centroid lies beyond the segment endpoint.
                endpoints=np.array([[float(row[f"ra{i}_deg"]),float(row[f"dec{i}_deg"])] for i in (0,1)])
                center=np.array([float(row["ra_deg"]),float(row["dec_deg"])])
                scale=np.array([np.cos(np.radians(center[1])),1.])*3600
                direction=(endpoints[1]-endpoints[0])*scale
                length=np.linalg.norm(direction)
                direction/=length
                sidevec=np.array([-direction[1],direction[0]])
                t=(2*int(row["axis_host_end"])-1)*(length/2+float(row["axis_gap_arcsec"]))
                hostsky=center+(direction*t+sidevec*float(row["axis_offset_arcsec"]))/scale
                host=wcs.world_to_pixel_values(*hostsky)
                host_label=f'axis gap {float(row["axis_gap_arcsec"]):.2f}"' 
            if vis[0].header["CD1_1"]>=0 or vis[0].header["CD2_2"]<=0:
                raise ValueError("unexpected native orientation")
            e0=np.array([float(row["x0"]),float(row["y0"])])
            e1=np.array([float(row["x1"]),float(row["y1"])])
            cen=(e0+e1)/2
            h,w=im.shape
            zoom=min(w,h,max(100,int(np.ceil(np.linalg.norm(e0-e1)))+34))
            left=int(np.clip(cen[0]-zoom/2,0,w-zoom))
            bot=int(np.clip(cen[1]-zoom/2,0,h-zoom))
            for j,bounds in enumerate([(0,0,w,h),(left,bot,zoom,zoom)]):
                sub=render(im,valid,bounds,side,(e0,e1,host),flip=True)
                page.paste(sub,(xx+4+j*(side+2),yy+20))
        draw.text((xx+4,yy+2),row["id"],fill="black")
        label=f'length {float(row["length_arcsec"]):.1f}"   {host_label}   fields {w*.1:.1f}" / {zoom*.1:.1f}"'
        draw.text((xx+4,yy+side+24),label,fill="black")
        if index==cols*rows_per_page-1 or i==len(targets)-1:
            path=a.out/f"vis_review_{i//(cols*rows_per_page)+1:02d}.jpg"
            used_rows=index//cols+1
            page.crop((0,0,page.width,26+ch*used_rows)).save(path,quality=88)
            paths.append(path.name)
    print(f"{len(targets)} segments on {len(paths)} pages")


if __name__=="__main__":
    main()
