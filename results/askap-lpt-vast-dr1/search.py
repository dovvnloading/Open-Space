#!/usr/bin/env python3
"""VAST Extragalactic DR1 LPT screen; catalogue epochs are NOT pulse periods."""
import argparse,hashlib,io,json,math,re,sys,traceback
from pathlib import Path
from urllib.parse import quote
import numpy as np
import pandas as pd
import requests

DOI="10.25919/nh9d-t846"
DAP="https://data.csiro.au/dap/ws/v2"
FILES=("VAST_Extragal_DR1_image_df.csv","VAST_Extragal_DR1_source_df.parquet","VAST_Extragal_DR1_measurement_df.parquet")
KNOWN=[
 ("GLEAM-X J162759.5-523504.3",246.997917,-52.584528,"precise"),
 ("GPM J1839-10",279.758333,-10.530381,"precise"),
 ("ASKAP J1832-0911",278.201912,-9.187583,"precise"),
 ("ASKAP J1935+2148",293.771358,21.811402,"precise"),
 ("ILT J1101+5521",165.460417,55.355444,"precise"),
 ("CHIME J0630+25",97.5,25.0,"rounded-name")]
CUTS={"release":DOI,"compactness":[0.75,1.4],"max_snr_min":8,"max_peak_mjy_min":5,
 "measurements_min":5,"selavy_min":2,"siblings_max":0,"positive_peak_ratio_min":5,
 "alt_eta_min":100,"alt_v_min":.8,"max_rank":60,"gaia_checks_max":12,
 "gaia_radius_arcsec":5,"atnf_radius_arcsec":30}
def digest(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for x in iter(lambda:f.read(8*1024*1024),b""):h.update(x)
 return h.hexdigest()
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(x,indent=2,sort_keys=True,default=str)+"\n")
def dist(a,b,c,d):
 a,b,c,d=[math.radians(float(z)) for z in [a,b,c,d]]
 z=math.sin((d-b)/2)**2+math.cos(b)*math.cos(d)*math.sin((c-a)/2)**2
 return math.degrees(2*math.asin(math.sqrt(max(0,min(1,z)))))*3600
def num(df,key,default=np.nan):
 return pd.to_numeric(df[key],errors="coerce") if key in df else pd.Series(default,index=df.index)
def screen(src):
 compact=num(src,"avg_compactness");n=num(src,"n_measurements",0);det=num(src,"n_selavy",0)
 sib=num(src,"n_siblings",0);snr=num(src,"max_snr",0)
 mx=num(src,"max_flux_peak");mn=num(src,"min_flux_peak");eta=num(src,"eta_peak",0);v=num(src,"v_peak",0)
 ratio=mx/np.maximum(mn,.5)
 extreme=((mn>0)&(ratio>=5))|((eta>=100)&(v>=.8))
 mask=compact.between(.75,1.4)&(n>=5)&(det>=2)&(sib<=0)&(snr>=8)&(mx>=5)&extreme
 result=src.loc[mask].copy()
 result["screen_ratio"]=ratio.loc[result.index]
 result["screen_score"]=(np.log1p(eta.loc[result.index].clip(lower=0))*(1+v.loc[result.index].clip(lower=0))
  *np.log1p(snr.loc[result.index].clip(lower=0))*np.log1p(mx.loc[result.index].clip(lower=0)))
 return result.sort_values(["screen_score","source"],ascending=[False,True])
def getjson(s,url,receipts,params=None):
 r=s.get(url,params=params,timeout=(15,70))
 receipts.append({"url":r.url,"status":r.status_code,"size":len(r.content),
                  "sha256":hashlib.sha256(r.content).hexdigest(),"type":"metadata"})
 r.raise_for_status();return r.json()
def acquire(s,cache,receipts):
 """Resolve the DOI to specific file objects, preserving all attempted routes."""
 cache.mkdir(parents=True,exist_ok=True)
 meta=None
 for ident in [DOI,"csiro:73391"]:
  try:meta=getjson(s,f"{DAP}/collections/{ident}",receipts);break
  except Exception:pass
 if meta is None:raise RuntimeError("CSIRO DAP collection metadata inaccessible")
 fed=meta.get("id",{}).get("identifier")
 if not fed:raise RuntimeError("CSIRO DAP collection identifier absent")
 write(cache/"collection_metadata.json",meta)
 receipts.append({"type":"collection_shape","keys":list(meta),
   "file_related_keys":[k for k in meta if "file" in k.lower() or "folder" in k.lower()],
   "identifiers":[str(meta.get("id",{}))[:250]],"frozen_metadata_sha256":digest(cache/"collection_metadata.json")})
 ids=[fed,quote(fed,safe=""),DOI,quote(DOI,safe=""),"102.100.100/713166",quote("102.100.100/713166",safe="")]
 listing=None;good_id=None
 for identifier in dict.fromkeys(ids):
  for extra in ["files.json","files"]:
   try:
    candidate=getjson(s,f"{DAP}/collections/{identifier}/{extra}",receipts,
       {"folder":"/lightcurves","page":1,"sb":"filename","size":50})
    if candidate.get("file") or candidate.get("files"):
     listing=candidate;good_id=identifier;break
   except Exception:pass
  if listing is not None:break
 if listing is not None:
  entries=listing.get("file",listing.get("files",[]))
  if isinstance(entries,dict):entries=[entries]
  urls={}
  for e in entries:
   raw=str(e.get("filename","")).lstrip("/")
   if raw.startswith("lightcurves/") and raw.split("/")[-1] in FILES:
    urls[raw.split("/")[-1]]=f"{DAP}/collections/{good_id}/downloadfile?fileName={quote('/'+raw,safe='')}"
  if set(urls)!=set(FILES):
   raise RuntimeError(f"DAP file index missing {set(FILES)-set(urls)}; seen {[e.get('filename') for e in entries[:8]]}")
 else:
  # The pinned upstream example notebook gives precisely these published
  # lightcurves file names. Test the image CSV before reading any large blob.
  urls=None
  for identifier in dict.fromkeys(ids):
   first=f"{DAP}/collections/{identifier}/downloadfile?fileName={quote('/lightcurves/'+FILES[0],safe='')}"
   try:
    with s.get(first,stream=True,timeout=(10,30)) as r:
     head=next(r.iter_content(256),b"")
     receipts.append({"type":"direct_download_probe","url":first,"http":r.status_code,
       "content_type":r.headers.get("content-type"),"head_sha256":hashlib.sha256(head).hexdigest(),
       "head_hex":head[:14].hex()})
     if r.status_code==200 and b"DOCTYPE" not in head and b"<html" not in head.lower():
      urls={name:f"{DAP}/collections/{identifier}/downloadfile?fileName={quote('/lightcurves/'+name,safe='')}" for name in FILES}
      break
   except Exception as e:
    receipts.append({"type":"direct_download_probe","url":first,"error":str(e)[:150]})
  if urls is None:raise RuntimeError("Collection metadata accessible, but all DAP files.json and direct-download routes failed")
 paths={}
 for name in FILES:
  p=cache/name
  if not p.exists() or p.stat().st_size<1024:
   part=Path(str(p)+".part");count=0
   with s.get(urls[name],timeout=(20,180),stream=True) as r:
    code=r.status_code;ctype=r.headers.get("content-type","");r.raise_for_status()
    with part.open("wb") as f:
     for b in r.iter_content(1024*1024):
      if not b:continue
      count+=len(b)
      if count>1500000000:raise RuntimeError(f"{name}: over 1.5GB")
      f.write(b)
   if count<1024:raise RuntimeError(f"{name}: download too small {count}")
   part.replace(p)
   receipts.append({"type":"input","url":urls[name],"file":name,"http":code,
      "content_type":ctype,"bytes":count,"sha256":digest(p)})
  else:receipts.append({"type":"input","url":urls[name],"file":name,
                        "http":"cached","bytes":p.stat().st_size,"sha256":digest(p)})
  if name.endswith(".parquet"):
   with p.open("rb") as f:
    if f.read(4)!=b"PAR1":raise RuntimeError(f"{name}: not parquet")
  paths[name]=p
 return paths,fed
def footprints(img):
 frame=img.copy()
 for c in ["centre_ra","centre_dec","physical_bmaj","physical_bmin"]:frame[c]=num(frame,c)
 frame=frame.loc[frame.physical_bmaj.between(.1,15)&frame.physical_bmin.between(.1,15)]
 rows=[]
 for name,ra,dec,precision in KNOWN:
  da=np.abs((frame.centre_ra-ra+180)%360-180)*np.cos(np.radians((frame.centre_dec+dec)/2)).clip(lower=.05)
  db=np.abs(frame.centre_dec-dec)
  hit=frame.loc[(da<=frame.physical_bmaj/2)&(db<=frame.physical_bmin/2)]
  rows.append({"name":name,"ra":ra,"dec":dec,"precision":precision,
               "images":[str(q) for q in hit.image_id],"n_square_envelope_images":len(hit),
               "qualifier":"not measured flux/sensitivity"})
 return rows,len(frame)
def epochs(parquet,leaders,out):
 import pyarrow.parquet as pq
 pf=pq.ParquetFile(parquet)
 cols=[c for c in ["source","image_id","time","flux_peak","flux_peak_err","flux_int",
 "flux_int_err","forced","snr","flag_c4","has_siblings","ra","dec","compactness"] if c in pf.schema.names]
 targets={int(x["source"]) for x in leaders}
 frames=[]
 for batch in pf.iter_batches(batch_size=150000,columns=cols):
  t=batch.to_pandas();t=t.loc[t.source.isin(targets)]
  if len(t):frames.append(t)
 if not frames:return {"status":"no_ranked_measurements","rows":0}
 m=pd.concat(frames,ignore_index=True).sort_values(["source","time"])
 m.to_csv(out/"ranked_epoch_measurements.csv",index=False)
 checks=[]
 for source,g in m.groupby("source"):
  d=g.loc[~g.forced.fillna(True).astype(bool)] if "forced" in g else g
  good=d.loc[~d.flag_c4.fillna(False).astype(bool)] if "flag_c4" in d else d
  flux=num(good,"flux_peak").clip(lower=0);err=num(good,"flux_peak_err",0).fillna(0)
  eff=np.sqrt(err**2+(.07*flux)**2)
  checks.append({"source":str(source),"epochs":len(g),"selavy":len(d),"good_selavy":len(good),
   "above_8sigma_including_7pct_systematic":int(((flux/eff.replace(0,np.nan))>=8).sum()),
   "first":str(g.time.min()),"last":str(g.time.max())})
 write(out/"epoch_quality.json",checks)
 return {"status":"measured","rows":len(m),"sources":len(checks)}
def gaia(s,leaders,out):
 url="https://gea.esac.esa.int/tap-server/tap/sync"
 result={"endpoint":url,"query_radius_arcsec":5,"checked":[],"unqueried":[]}
 for i,item in enumerate(leaders):
  ident=str(item["source"])
  if i>=12:result["unqueried"].append(ident);continue
  ra,dec=float(item["wavg_ra"]),float(item["wavg_dec"])
  q=("SELECT TOP 20 source_id,ra,dec,parallax,phot_g_mean_mag FROM gaiadr3.gaia_source "
     f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},{5/3600:.9f}))")
  row={"source":ident,"state":"unavailable"}
  try:
   r=s.get(url,params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":q},timeout=(10,30))
   row.update({"status":r.status_code,"sha256":hashlib.sha256(r.content).hexdigest()});r.raise_for_status()
   table=pd.read_csv(io.StringIO(r.text),comment="#")
   if not {"ra","dec"}.issubset(table.columns):raise RuntimeError("Gaia CSV missing coords")
   row["matches"]=[{"gaia_id":str(x["source_id"]),
      "arcsec":round(dist(ra,dec,x["ra"],x["dec"]),2),"parallax_mas":None if pd.isna(x.parallax) else float(x.parallax)}
      for _,x in table.iterrows() if dist(ra,dec,x["ra"],x["dec"])<=5]
   row["state"]="queried"
  except Exception as e:row["error"]=str(e)[:180]
  result["checked"].append(row)
 write(out/"gaia_receipts.json",result);return result
def atnf(s,leaders,out):
 url="https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat.db"
 status={"url":url,"state":"unavailable","nearest":[]}
 try:
  r=s.get(url,timeout=(20,75));r.raise_for_status()
  status.update({"http":r.status_code,"bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest()})
  from astropy.coordinates import SkyCoord
  import astropy.units as u
  pulsars=[]
  for block in r.text.split("@"):
   props=dict(re.findall(r"(?m)^(PSRJ|RAJ|DECJ)\s+(\S+)",block))
   if set(("PSRJ","RAJ","DECJ"))<=set(props):
    try:
     c=SkyCoord(props["RAJ"],props["DECJ"],unit=(u.hourangle,u.deg))
     pulsars.append((props["PSRJ"],c.ra.deg,c.dec.deg))
    except Exception:pass
  if len(pulsars)<100:raise RuntimeError(f"Only {len(pulsars)} parsed pulsars")
  status["state"]="queried";status["n_pulsars"]=len(pulsars)
  for lead in leaders:
   d,n=min((dist(lead["wavg_ra"],lead["wavg_dec"],p[1],p[2]),p[0]) for p in pulsars)
   status["nearest"].append({"source":str(lead["source"]),"psr":n,
      "arcsec":round(d,2),"within_30arcsec":d<=30})
 except Exception as e:status["error"]=str(e)[:200]
 write(out/"atnf_receipts.json",status);return status
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
 p.add_argument("--cache",type=Path,default=Path("/tmp/vast-dr1"));args=p.parse_args()
 out=args.out;out.mkdir(parents=True,exist_ok=True)
 write(out/"frozen_protocol.json",CUTS)
 summary={"state":"running","claim":"NO VALIDATED NEW LONG-PERIOD TRANSIENT",
          "not_a_period_search":True,"protocol_sha256":digest(out/"frozen_protocol.json")}
 receipts=[];s=requests.Session();s.headers["User-Agent"]="Open-Space-VAST-DR1-LPT-research/1.0"
 try:
  files,fed=acquire(s,args.cache,receipts);summary["dap_collection"]=fed
  images=pd.read_csv(files[FILES[0]])
  if "image_id" not in images and "Unnamed: 0" in images:
   images.rename(columns={"Unnamed: 0":"image_id"},inplace=True)
  if "image_id" not in images:raise RuntimeError(f"No image_id in {list(images.columns)}")
  images.to_json(out/"frozen_images.jsonl",orient="records",lines=True,date_format="iso")
  summary.update({"frozen_image_sha256":digest(out/"frozen_images.jsonl"),
                  "frozen_images":len(images),"distinct_images":int(images.image_id.nunique())})
  controls,n=footprints(images)
  write(out/"known_lpt_footprints.json",{"controls":controls,"geometry_eligible_images":n})
  summary["controls_in_approximate_footprints"]=sum(x["n_square_envelope_images"]>0 for x in controls)
  src=pd.read_parquet(files[FILES[1]])
  if "source" not in src and src.index.name=="source":src=src.reset_index()
  summary["source_rows"]=len(src)
  rank=screen(src);summary["screen_all"]=len(rank);rank=rank.head(60)
  summary["ranked_rows"]=len(rank);rank.to_csv(out/"ranked_sources.csv",index=False)
  import pyarrow.parquet as pq
  summary["measurement_rows"]=pq.ParquetFile(files[FILES[2]]).metadata.num_rows
  leaders=rank.to_dict("records")
  summary["epoch_review"]=epochs(files[FILES[2]],leaders,out) if leaders else {"status":"no_leads","rows":0}
  summary["gaia"]=gaia(s,leaders,out) if leaders else {"status":"no_leads"}
  summary["atnf"]=atnf(s,leaders,out) if leaders else {"status":"no_leads"}
  summary["contaminant_prior"]={"published_2026_DR1_variability_classes":{"AGN":44,"stars":40,"pulsars":27,"other":6},
    "AGN_and_stars_fraction_of_published_117":84/117,
    "false_positive_rate":"NOT CALIBRATED; published source classes are not an independent validation set"}
  summary["pending"]=["VAST native FITS cutout inspection and field-specific noise","RACS low/mid/high I/V/Q/U photometry",
   "Legacy Surveys/unWISE counterpart exclusion","10-second or snapshot periodicity (not in released DR1 tables)",
   "eRASS1/Swift X-ray, SIMBAD/literature checks","empirical AGN scintillation/flare-star false-positive calibration"]
  summary["state"]="catalogue_screen_only_no_periodicity"
 except Exception as e:
  summary.update({"state":"incomplete_access_or_analysis","error":f"{type(e).__name__}: {e}",
    "traceback":traceback.format_exc(limit=4)})
 finally:
  write(out/"network_receipts.json",receipts);write(out/"summary.json",summary)
  print(json.dumps({k:v for k,v in summary.items() if k not in ["gaia","atnf","pending","traceback"]},indent=2,default=str))
 return 1 if summary["state"]=="incomplete_access_or_analysis" else 0
if __name__=="__main__":sys.exit(main())
