"""Blended forced photometry at five additional catalogue positions in IRSA Q1 maps."""
import argparse
import csv
import hashlib
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from cutouts import image_extensions
from measure import HERE,blank_measurement,ab,stamp,table

IDS=("TILE102021981_-584133901477397990","TILE102160059_2706053122676895502",
     "TILE102020063_-624125255493125274","TILE102159491_2752873964662719157",
     "TILE102019126_-587713107502760148")


def pix(h,ra,dec):
    ra,dec,ra0,d0=np.deg2rad([ra,dec,float(h["CRVAL1"]),float(h["CRVAL2"])])
    den=np.sin(dec)*np.sin(d0)+np.cos(dec)*np.cos(d0)*np.cos(ra-ra0)
    tangent=np.rad2deg([np.cos(dec)*np.sin(ra-ra0)/den,
                        (np.sin(dec)*np.cos(d0)-np.cos(dec)*np.sin(d0)*np.cos(ra-ra0))/den])
    cd=[[float(h["CD1_1"]),float(h["CD1_2"])],[float(h["CD2_1"]),float(h["CD2_2"])]]
    return np.linalg.solve(cd,tangent)+np.array([float(h["CRPIX1"]),float(h["CRPIX2"])])-1


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("agile_tsv",type=Path);args=p.parse_args()
    lines=args.agile_tsv.read_text().splitlines();index=next(i for i,x in enumerate(lines) if x.startswith("recno\t"))
    catalog={r["cutoutName"].strip():r for r in csv.DictReader(lines[index:],delimiter="\t") if r.get("cutoutName")}
    def analyze(obj):
        r=catalog[obj];ra=float(r["RAJ2000"]);dec=float(r["DEJ2000"])
        tile=obj.split("_")[0][4:]; phot,manifest,draw=[],[],[]
        for b in ("VIS","H"):
            sec,key=("VIS","VIS") if b=="VIS" else ("NISP","NIR-H")
            directory=f"https://irsa.ipac.caltech.edu/ibe/data/euclid/q1/MER/{tile}/{sec}/"
            listing=urlopen(directory,timeout=60).read().decode()
            arrays=[]
            for product,prefix in [("flux",f"EUC_MER_BGSUB-MOSAIC-{key}_TILE"),("rms",f"EUC_MER_MOSAIC-{key}-RMS_TILE")]:
                filenames=re.findall(r'href="('+prefix+r'[^"]+\.fits)"',listing)
                if len(filenames)!=1:raise ValueError((obj,b,filenames))
                url=directory+filenames[0]+f"?center={ra},{dec}deg&size=30arcsec&gzip=false"
                data=urlopen(url,timeout=90).read()
                hdr,im=next(iter(image_extensions(data).values()));arrays.append(im)
                if product=="flux":h=hdr
                manifest.append(dict(id_str=obj,band=b,product=product,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),url=url))
            im,rms=arrays;x,y=pix(h,ra,dec);valid=np.isfinite(im)&np.isfinite(rms)&(rms>0)&(rms<1e15)
            dy,dx,bg,sig,n=blank_measurement(im,rms,valid,y,x,3)
            ap=(round(y)+dy,round(x)+dx);counts=im[ap].sum()-bg
            sigma=max(sig,np.sqrt(np.sum(np.square(rms[ap].astype(float)))))
            if not valid[ap].all():counts=sigma=np.nan
            phot.append(dict(id_str=obj,ra_catalogue=ra,dec_catalogue=dec,grade=r["Grade"],band=b,
                  catalogue_core_counts=counts,noise=sigma,blended_snr=counts/sigma,
                  blended_mag=ab(float(h["MAGZERO"]),counts),noise_only_5sigma_mag=ab(float(h["MAGZERO"]),5*sigma)))
            draw.append((b+" observed",im[round(y)-35:round(y)+36,round(x)-35:round(x)+36]))
        stamp(draw,obj,HERE)
        return phot,manifest
    phot,manifest=[],[]
    with ThreadPoolExecutor(5) as pool:
        for f,m in pool.map(analyze,IDS):phot.extend(f);manifest.extend(m)
    table(HERE/"agile_followup.csv",phot);table(HERE/"agile_followup_inputs.csv",manifest)


if __name__=="__main__":main()
