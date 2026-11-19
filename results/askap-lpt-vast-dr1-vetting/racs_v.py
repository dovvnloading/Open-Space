#!/usr/bin/env python3
"""Bounded probe of the public RACS-low2 I/V catalogue release, not an image-based polarimetry validation."""
import csv,hashlib,io,json,math,re
from urllib.parse import quote
import requests
import numpy as np
from astropy.table import Table
BASE="https://data.csiro.au/dap/ws/v2"
DOI="10.25919/m8t2-b486"
def _json(url,params=None):
 r=requests.get(url,params=params,timeout=(12,38),headers={"User-Agent":"Open-Space-RACS-V-followup/1.0"})
 r.raise_for_status()
 return r.json()
def probe(positions):
 """Return positive 221-source catalogue matches if accessible; absence NEVER means Stokes V nondetection."""
 out={"release_doi":DOI,"status":"unavailable","images_not_retrieved":True,"not_a_stokes_v_upper_limit":True,
      "matches":{},"checked":len(positions),"directory_listing":[]}
 try:
  meta=_json(f"{BASE}/collections/{DOI}")
  fid=meta.get("id",{}).get("identifier","csiro:74963")
  out["fedora_id"]=fid
  filelist=[]
  # A small, bounded inventory: do not attempt massive image release traversal.
  for folder in ("/","/catalogues","/catalogs","/catalogue","/data"):
   try:
    obj=_json(f"{BASE}/collections/{fid}/files",
              params={"folder":folder,"page":1,"size":100,"sb":"filename"})
    files=obj.get("file",obj.get("files",[]))
    if isinstance(files,dict):files=[files]
    for f in files[:100]:
     name=f.get("filename",f.get("name",""))
     if name and name not in [x.get("filename") for x in filelist]:
      filelist.append(f)
    out["directory_listing"].append({"folder":folder,"status":"ok","count_page":len(files)})
   except Exception as ex:
    out["directory_listing"].append({"folder":folder,"status":"unavailable","error":str(ex)[:110]})
  paths=[]
  for f in filelist:
   name=str(f.get("filename",f.get("name","")))
   suffix=name.lower().rsplit(".",1)[-1]
   if suffix not in ("csv","fits","fit","vot","xml","tsv","txt"):continue
   if "cat" not in name.lower() and "source" not in name.lower():continue
   if re.search(r"(stokes.?v|circular|racs[_-]?low2.?v|catalogue.?v)",name,re.I):
    paths.append(name)
  out["v_catalogue_candidates"]=paths[:12]
  for path in paths[:4]:
   url=f"{BASE}/collections/{fid}/downloadfile?fileName={quote('/'+path.lstrip('/'),safe='')}"
   try:
    r=requests.get(url,timeout=(12,80),stream=True,headers={"User-Agent":"Open-Space-RACS-V-followup/1.0"})
    r.raise_for_status(); b=io.BytesIO()
    for chunk in r.iter_content(1024*1024):
     b.write(chunk)
     if b.tell()>40_000_000:raise ValueError("catalogue exceeds 40MB cap")
    data=b.getvalue()
    if not data:raise ValueError("empty catalogue file")
    name=path.lower()
    if name.endswith(".csv"): t=Table.read(io.BytesIO(data),format="ascii.csv")
    elif name.endswith(".tsv"):t=Table.read(io.BytesIO(data),format="ascii.tab")
    elif name.endswith((".fits",".fit")):t=Table.read(io.BytesIO(data),format="fits")
    elif name.endswith((".vot",".xml")):t=Table.read(io.BytesIO(data),format="votable")
    else:t=Table.read(io.BytesIO(data),format="ascii")
    cols={str(x).lower():str(x) for x in t.colnames}
    a,b=next(((cols[x],cols[y]) for x,y in
          (("ra","dec"),("raj2000","dej2000"),("ra_deg","dec_deg"),("ra_icrs","de_icrs"))
          if x in cols and y in cols),(None,None))
    if not a:raise ValueError(f"no coordinate columns in {list(cols)[:12]}")
    if len(t)>1_000_000:raise ValueError("not small thresholded V catalogue")
    for sid,(ra,dec) in positions.items():
     cand=[]
     for row in t:
      try:
       dr=abs(float(row[a])-ra)*math.cos(math.radians(dec))
       dd=float(row[b])-dec
       sep=3600*math.hypot(dr,dd)
       if sep<20:cand.append((sep,row))
      except (TypeError,ValueError):continue
     if cand:
      sep,row=min(cand,key=lambda x:x[0])
      out["matches"][sid]={"arcsec":round(sep,2),"record":{x:str(row[x]) for x in t.colnames[:20]}}
    out.update({"status":"selected_catalogue_crossmatched","selected_catalogue":path,
       "bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),"rows":len(t),
       "matched_sources":len(out["matches"]),
       "interpretation":"positive matches only; no match is NOT a Stokes V nondetection; measure image V and leakage to validate"})
    break
   except Exception as ex:
    out.setdefault("catalogue_errors",[]).append({"path":path,"error":f"{type(ex).__name__}: {str(ex)[:160]}"})
  if not paths:out["status"]="no_small_V_catalogue_in_bounded_DAP_directory_listing"
 except Exception as ex:
  out["error"]=f"{type(ex).__name__}: {str(ex)[:200]}"
 return out
