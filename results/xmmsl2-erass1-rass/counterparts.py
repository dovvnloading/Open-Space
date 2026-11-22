"""Read the published Salvato et al. eRASS1 counterpart associations for the leads."""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import tarfile
import tempfile
from html.parser import HTMLParser
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from scipy.spatial import cKDTree

ROOT = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/SalvatoM_DR1/"
RELEASES = ("LS10","GDR3","CW2020")


class Text(HTMLParser):
    def __init__(self):
        super().__init__()
        self.text=[]
    def handle_data(self,value):
        if value.strip():
            self.text.append(value.strip())


def norm(x):
    return re.sub(r"[^a-z0-9]","",str(x).lower())


def as_text(x):
    if isinstance(x,bytes):
        return x.decode("utf-8","replace").strip().strip("\0")
    return str(x).strip().strip("\0")


def values(hdu,col):
    return np.asarray(hdu.data[col]).ravel()


def xyz(ra,dec):
    a,d=np.radians(ra),np.radians(dec)
    return np.column_stack((np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)))


def matched_rows(hdu, leads, note):
    cols=hdu.columns.names
    lookup={norm(k):k for k in cols}
    def first(names):
        return next((lookup[x] for x in names if x in lookup),None)
    # Exact eROSITA source-name join is preferred to a positional join.
    name=first(["eroiauname","iauname","erass1iauname","eroname","erass1name","erassname","erass1"])
    out={}
    if name:
        wanted={as_text(x).replace(" ",""):x for x in leads.ero_name}
        for i,s in enumerate(values(hdu,name)):
            key=as_text(s).replace(" ","")
            if key in wanted:
                out.setdefault(wanted[key],[]).append(i)
        if out:
            note.update(key=name,method="exact released eROSITA designation")
            return out
    # Require both coordinates of the eROSITA position. Do not substitute
    # an arbitrary nearby optical object for a released X-ray association.
    for rx,dx in [("erora","erodec"),("raero","decero"),("raerass1","decerass1"),
                  ("erass1ra","erass1dec"),("rax","decx"),("xrayra","xraydec"),("ra","dec")]:
        if rx not in lookup or dx not in lookup:
            continue
        x=values(hdu,lookup[rx]).astype(float)
        d=values(hdu,lookup[dx]).astype(float)
        good=np.isfinite(x)&np.isfinite(d)&(np.abs(d)<=90)&(x>=0)&(x<360)
        if not good.any():
            continue
        valid=np.flatnonzero(good)
        tree=cKDTree(xyz(x[good],d[good]))
        # A tenth of an arcsecond tests identity with a catalogue coordinate,
        # well inside the random error circle of a counterpart association.
        ix=tree.query_ball_point(xyz(leads.ero_ra,leads.ero_dec),2*np.sin(np.radians(.1/3600)/2))
        mapping={n:valid[v].tolist() for n,v in zip(leads.ero_name,ix) if len(v)}
        if mapping:
            note.update(key=[lookup[rx],lookup[dx]],method="agreement with released eROSITA position within 0.1 arcsec")
            return mapping
    raise ValueError("No exact eROSITA identifier or catalogue-coordinate join; inspect column definitions")


def read_release(tag, leads, out):
    base="Salvato_etal2025_DR1_"+tag
    url=ROOT+base+".colfits.tgz"
    note={"catalog":tag,"url":url}
    rows=[]
    try:
        response=requests.get(ROOT+base+".html",timeout=(15,120))
        response.raise_for_status()
        p=Text()
        p.feed(response.text)
        note["documentation_url"]=response.url
        note["documentation_sha256"]=hashlib.sha256(response.content).hexdigest()
        (out/("model_"+tag+".txt")).write_text("\n".join(p.text))
        with tempfile.TemporaryDirectory() as root:
            package=Path(root)/"release.tgz"
            with requests.get(url,stream=True,timeout=(15,180)) as r, package.open("wb") as dest:
                r.raise_for_status()
                for block in r.iter_content(2**20):
                    dest.write(block)
            with package.open("rb") as dest:
                note["sha256"]=hashlib.file_digest(dest,"sha256").hexdigest()
            note["bytes"]=package.stat().st_size
            with tarfile.open(package,"r|gz") as tf:
                for m in tf:
                    if m.isfile() and m.name.lower().endswith((".fits",".fit",".fits.gz")):
                        path=Path(root)/("table.fits.gz" if m.name.endswith(".gz") else "table.fits")
                        with path.open("wb") as dest:
                            shutil.copyfileobj(tf.extractfile(m),dest)
                        break
                else:
                    raise ValueError("No FITS table in release")
            with fits.open(path,memmap=True,character_as_bytes=True) as hdus:
                hdu=next(h for h in hdus if isinstance(h,fits.BinTableHDU))
                cols=hdu.columns.names
                note["fields"]=[(c.name,str(c.unit or ""),str(c.format),str(c.dim or "")) for c in hdu.columns]
                mapping=matched_rows(hdu,leads,note)
                # Retain the published fields without imposing a generic cut on
                # probabilities that may use different definitions across releases.
                for names,idxs in mapping.items():
                    for i in idxs:
                        record=dict(ero_name=names,catalog=tag,release_row=i+1)
                        for col in cols:
                            field=values(hdu,col)
                            if i<len(field):
                                record["pub_"+col]=as_text(field[i])
                        rows.append(record)
        note.update(status="ok",matched_ero_names=len(mapping),matched_rows=len(rows))
        print(tag, "join:", note["key"], "fields:", [x[0] for x in note["fields"]],flush=True)
    except Exception as e:
        note["error"]=repr(e)
    print(tag,note.get("matched_rows"),note.get("error",""),flush=True)
    return rows,note


def hard_catalog(inputs,dest):
    url=ROOT.rsplit("SalvatoM_DR1/",1)[0]+"MerloniA_DR1/eRASS1_Hard.tar.gz"
    result=[]
    note={"url":url,"band_keV":[2.3,5.0]}
    try:
        from astropy.coordinates import SkyCoord, search_around_sky
        import astropy.units as u
        response=requests.get(url,timeout=(15,120))
        response.raise_for_status()
        note["sha256"]=hashlib.sha256(response.content).hexdigest()
        with tarfile.open(fileobj=__import__("io").BytesIO(response.content)) as tar:
            member=next(m for m in tar if m.isfile() and m.name.lower().endswith(".fits"))
            with fits.open(tar.extractfile(member),memmap=False) as hdus:
                h=next(h for h in hdus if isinstance(h,fits.BinTableHDU))
                cols=set(h.columns.names)
                pos=SkyCoord(h.data["RA"]*u.deg,h.data["DEC"]*u.deg)
                q=SkyCoord(inputs.ra.to_numpy()*u.deg,inputs.dec.to_numpy()*u.deg)
                a,b,sep,_=search_around_sky(q,pos,35*u.arcsec)
                note["catalogue_rows"]=len(pos)
                for i,j,d in zip(a,b,sep.arcsec):
                    row={"detection_uid":inputs.iloc[i].detection_uid,"separation_arcsec":d}
                    for c in ("IAUNAME","RA","DEC","RADEC_ERR","DET_LIKE_0","DET_LIKE_3","ML_FLUX_3","ML_FLUX_ERR_3","MJD"):
                        if c in cols:
                            row[c]=as_text(h.data[c][j])
                    result.append(row)
        note["neighbors_within_35arcsec"]=len(result)
    except Exception as e:
        note["error"]=repr(e)
    pd.DataFrame(result,columns=["detection_uid","separation_arcsec","IAUNAME","RA","DEC","RADEC_ERR","DET_LIKE_0","DET_LIKE_3","ML_FLUX_3","ML_FLUX_ERR_3","MJD"]).to_csv(dest/"erass_hard.csv",index=False)
    return note


def main(out,source_file):
    inputs=pd.read_csv(out/source_file,low_memory=False)
    if "id" not in inputs:
        inputs["id"]=inputs.detection_uid
    sources=inputs.dropna(subset=["ero_name","ero_ra","ero_dec"]).drop_duplicates("ero_name")
    dest=out/"published_hosts"
    dest.mkdir(exist_ok=True)
    audit=[hard_catalog(inputs,dest)]
    for tag in RELEASES:
        rows,rec=read_release(tag,sources,dest) if len(sources) else ([],{"catalog":tag,"status":"no eROSITA counterparts in lead set"})
        audit.append(rec)
        pd.DataFrame(rows).to_csv(dest/(tag+".csv.gz"),index=False)
    (dest/"manifest.json").write_text(json.dumps(audit,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    p.add_argument("--source-file",default="input.csv")
    a=p.parse_args()
    main(a.out,a.source_file)
