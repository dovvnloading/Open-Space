#!/usr/bin/env python3
"""Metadata-only S3 discovery when IRSA SIA is returning HTTP 502.
This pre-pixel step may not be interpreted as a nondetection.
"""
from pathlib import Path
from collections import Counter
from urllib.parse import quote
from xml.etree import ElementTree as ET
import csv,io,json,math,hashlib
import requests
import s3fs
import numpy as np
import pyarrow.parquet as pq
import pyarrow.fs as pafs
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
P=Path(__file__).resolve().parent
O=P/"output"
BUCKET="nasa-irsa-euclid-dr1"
BASE="https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/"
PREFIX="contributed/dr1/hats/reformatted_catalogs/mer_catalogue_deep_survey/euclid_dr1_reformatted_catalogs_mer_catalogue_deep_survey-hats/"
FIELDS=json.loads((P/"freeze.json").read_text())["fields"]
def sha(s):return hashlib.sha256(s).hexdigest()
def req(url):
    r=requests.get(url,timeout=(12,80));r.raise_for_status();return r.content
def list_prefix(prefix):
    raw=req(BASE+"?list-type=2&prefix="+quote(prefix,safe="/"))
    root=ET.fromstring(raw)
    keys=[x.text for x in root.findall(".//{*}Contents/{*}Key") if x.text]
    return keys,{"query_prefix":prefix,"listing_sha256":sha(raw),"key_count_first_page":len(keys)}
def one(field,index,fs,arrow):
    f=field["name"];ra,dec=field["ra_deg"],field["dec_deg"]
    out={"field":f,"center":[ra,dec],"status":"PENDING","hats_catalogue":PREFIX,"tiles":[]}
    # A HATS leaf is preselected from its spatial partition index, not from
    # a transient-like pixel. Read only 3 nearby leaves when necessary.
    leaves=[]
    for x in index:
        o,p=int(x["Norder"]),int(x["Npix"])
        h=HEALPix(nside=2**o,order="nested",frame="icrs")
        same=int(h.lonlat_to_healpix(ra*u.deg,dec*u.deg))==p
        if same:leaves.append(x)
    if not leaves:
        out["error"]="NO_HATS_LEAF_CONTAINING_FIELD_CENTER"
        return out
    out["selected_hats_leaves"]=[{k:v for k,v in e.items() if k in ("Norder","Npix","Dir","num_rows")} for e in leaves]
    hits=[]
    for e in leaves[:3]:
        folder=PREFIX+f"dataset/Norder={e['Norder']}/Dir={e['Dir']}/Npix={e['Npix']}/"
        files,lr=list_prefix(folder)
        out.setdefault("hats_file_receipts",[]).append(lr)
        names=[k for k in files if k.endswith(".parquet")]
        if len(names)!=1:raise ValueError("expected one Parquet file per HATS partition "+str(names))
        with arrow.open_input_file(BUCKET+"/"+names[0]) as h:
            pf=pq.ParquetFile(h)
            cols=pf.schema_arrow.names
            idcol=next((c for c in ("tileid","tile_id","TILEID") if c in cols),None)
            if not idcol:raise ValueError("No tile identifier among columns: "+str(cols[:70]))
            out.setdefault("hats_columns",list(cols[:50]))
            for j in range(pf.num_row_groups):
                tab=pf.read_row_group(j,columns=["ra","dec",idcol])
                rr=np.asarray(tab["ra"].to_numpy(zero_copy_only=False),float)
                dd=np.asarray(tab["dec"].to_numpy(zero_copy_only=False),float)
                dis=np.hypot((rr-ra)*math.cos(math.radians(dec)),dd-dec)
                ix=np.where(np.isfinite(dis)&(dis<.05))[0]
                for k in ix:
                    hits.append((float(dis[k]),str(int(tab[idcol][int(k)].as_py()))))
    out["hats_sources_within_3arcmin"]=len(hits)
    counts=Counter(t for _,t in hits)
    out["tile_counts"]=dict(counts.most_common(8))
    tileids=[x for x,_ in counts.most_common(2)]
    if f=="north" and "102158276" not in tileids:tileids=["102158276"]+tileids[:1]
    if not tileids:out["status"]="NO_CATALOG_TILE_NEAR_FROZEN_CENTER";return out
    for tile in tileids[:2]:
        part={"tileid":tile,"prefixes":[],"samples":[]}
        for kind in ("VIS","NISP"):
            prefix=f"dr1/MER/{tile}/{kind}/"
            try:
                keys,lr=list_prefix(prefix)
                part["prefixes"].append(lr)
                categorized={}
                for key in keys:
                    name=key.rsplit("/",1)[-1]
                    cat=next((name2 for name2,token in (("science","BGSUB-MOSAIC"),("rms","-RMS_"),
                        ("flag","-FLAG_"),("layer","DETECTOR-LAYERING")) if token in name),None)
                    if cat:
                        band="VIS" if kind=="VIS" else next((b for b in ("Y","J","H") if ("NIR-"+b+"_TILE") in name or ("NIR-"+b+"-") in name),None)
                        if band:categorized.setdefault(band+"_"+cat,[]).append(key)
                part.setdefault("counts",{}).update({k:len(v) for k,v in categorized.items()})
                part.setdefault("files",{}).update({k:sorted(v)[:8] for k,v in categorized.items()})
                for ident,keys in categorized.items():
                    if not ident.endswith("_science"):continue
                    for key in sorted(keys)[:3]:
                        try:
                            src="s3://"+BUCKET+"/"+key
                            with fits.open(src,use_fsspec=True,
                                fsspec_kwargs={"anon":True,"client_kwargs":{"region_name":"us-east-1"}}) as hd:
                                hdr=hd[0].header
                                wanted=["OBS_ID","OBSID","OBSERVAT","VISITID","OBSID","MJD-OBS","DATE-OBS",
                                    "DATE","CRVAL1","CRVAL2","NAXIS1","NAXIS2","FILENAME","EQUINOX","FILTER",
                                    "TELESCOP","PRODTYPE","BUNIT"]
                                cards={k:str(hdr.get(k)) for k in wanted if k in hdr}
                                part["samples"].append({"key":key,"hdu_count":len(hd),"header":cards})
                        except Exception as e:part["samples"].append({"key":key,"error":type(e).__name__+": "+str(e)[:240]})
            except Exception as e:
                part.setdefault("errors",{})[kind]=type(e).__name__+": "+str(e)[:240]
        out["tiles"].append(part)
    out["status"]="S3_DISCOVERY_PROBE_COMPLETE"
    return out
def main():
    O.mkdir(exist_ok=True)
    out={"status":"S3_PROBE_IN_PROGRESS","fields":[]}
    try:
        blob=req(BASE+PREFIX+"partition_info.csv")
        out["hats_index_sha256"]=sha(blob);out["hats_index_bytes"]=len(blob)
        index=list(csv.DictReader(io.StringIO(blob.decode("utf-8-sig"))))
        out["hats_partitions"]=len(index)
        fs=s3fs.S3FileSystem(anon=True,client_kwargs={"region_name":"us-east-1"})
        arrow=pafs.S3FileSystem(anonymous=True,region="us-east-1")
        for field in FIELDS:
            try:r=one(field,index,fs,arrow)
            except Exception as e:r={"field":field["name"],"status":"S3_PROBE_ERROR","error":type(e).__name__+": "+str(e)[:480]}
            out["fields"].append(r)
            (O/"s3_probe.json").write_text(json.dumps(out,indent=2)+"\n")
            print("S3_FIELD",r["field"],r["status"],[x["tileid"] for x in r.get("tiles",[])],flush=True)
    except Exception as e:
        out.update(status="S3_INDEX_ERROR",error=type(e).__name__+": "+str(e)[:400])
    if out["status"]!="S3_INDEX_ERROR":out["status"]="S3_PROBE_DONE"
    (O/"s3_probe.json").write_text(json.dumps(out,indent=2)+"\n")
    print("S3_PROBE",out["status"],flush=True)
if __name__=="__main__":main()
