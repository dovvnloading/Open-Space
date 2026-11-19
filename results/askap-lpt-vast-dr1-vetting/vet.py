#!/usr/bin/env python3
"""Native-image and counterpart follow-up of frozen VAST DR1 screen (PR 200)."""
import csv,hashlib,io,json,math,os,pickletools,sys,traceback
from pathlib import Path
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BASE=Path("results/askap-lpt-vast-dr1/output")
OUT=Path("results/askap-lpt-vast-dr1-vetting/output")
INDEX="https://raw.githubusercontent.com/IrisdeRuiter/VAST_Extragalactic_DR1/main/download_tools/cutouts_dict.pkl"
BLOB="b0dc867d6baf4a967d84fba781d2b27ed96bbe67"
CAT={
 "racs_low_i":("J/other/PASA/38.58",7),
 "racs_mid_i":("J/other/PASA/41.3/gcompsm",7),
 "racs_high_i":("J/other/PASA/42.38/gcompsh",7),
 "gaia":("I/355/gaiadr3",3),
 "legacy_dr10_photoz_subset":("VII/300",3),
 "unwise":("II/363/unwise",3),
 "panstarrs":("II/349/ps1",3),
 "erass1":("J/A+A/682/A34/erass1-m",12),
 "swift":("IX/58/2sxps",12)}
FREQ={"racs_low_i":887.5,"racs_mid_i":1367.5,"racs_high_i":1655.5}
def dump(name,x):
 OUT.mkdir(exist_ok=True,parents=True)
 (OUT/name).write_text(json.dumps(x,indent=2,sort_keys=True,default=str,allow_nan=False)+"\n")
def flag(x):return str(x).lower() in ("true","t","1")
def val(x):
 try:
  f=float(x)
  return f if math.isfinite(f) else None
 except (TypeError,ValueError):return None
def angular(a,b,c,d):
 a,b,c,d=map(math.radians,[a,b,c,d])
 h=math.sin((d-b)/2)**2+math.cos(b)*math.cos(d)*math.sin((c-a)/2)**2
 return 7200*math.degrees(math.asin(min(1,math.sqrt(max(0,h)))))
def shortlist():
 s=pd.read_csv(BASE/"ranked_sources.csv",dtype={"source":str})
 m=pd.read_csv(BASE/"ranked_epoch_measurements.csv",dtype={"source":str})
 assert len(s)==60 and len(m)==693 and m.source.nunique()==60
 on=[]
 for sid,g in m.groupby("source"):
  bright=g.loc[~g.forced.map(flag)&~g.flag_c4.map(flag)&(g.snr>=8)]
  faint=g.loc[g.forced.map(flag)&(g.snr.abs()<3)]
  if len(bright)>=2 and len(faint)>0:on.append(sid)
 assert len(on)==14,(len(on),on)
 ids=[x for x in s.source if x in on]+[x for x in s.source if x not in on]
 return s.set_index("source"),m,ids,on
def get_file(url,maxlen=40_000_000,seconds=110):
 r=requests.get(url,timeout=(12,seconds),stream=True,headers={"User-Agent":"Open-Space-LPT-followup/1.0"})
 r.raise_for_status();data=io.BytesIO()
 for b in r.iter_content(262144):
  data.write(b)
  if data.tell()>maxlen:raise RuntimeError("response byte limit")
 return data.getvalue(),r.url
def cutout_index(ids):
 temp=Path(os.getenv("RUNNER_TEMP","/tmp"))/"vast_dr1_public_index.pkl"
 receipt={"url":INDEX,"status":"failed","git_blob_verified":False,"matches":0}
 try:
  if not temp.exists():
   b,url=get_file(INDEX,maxlen=120_000_000,seconds=180)
   temp.write_bytes(b)
  b=temp.read_bytes()
  digest=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()
  receipt.update({"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),
                  "git_blob_verified":digest==BLOB})
  if digest!=BLOB:raise ValueError("upstream cutout dictionary changed since pin")
  found={};want=set(ids);pending=None
  with temp.open("rb") as f:
   # Never call pickle.load on a remotely downloaded pickle.
   for op,arg,offset in pickletools.genops(f):
    if isinstance(arg,str):
     if arg in want:pending=arg
     elif pending and arg.startswith("https://data.csiro.au/dap/ws/v2/collections/") and "/downloadfile?" in arg:
      found[pending]=arg;pending=None
  receipt.update({"status":"ok","matches":len(found),"missing":sorted(want-set(found))})
  return found,receipt
 except Exception as e:
  receipt["error"]=f"{type(e).__name__}: {str(e)[:160]}"
  return {},receipt
def pixels(arr,head,ra,dec):
 a=np.asarray(arr,dtype=float).squeeze()
 if a.ndim!=2 or min(a.shape)<12:return {"state":"invalid_shape"},a
 ny,nx=a.shape
 try:
  x,y=WCS(head).celestial.world_to_pixel(SkyCoord(ra*u.deg,dec*u.deg))
  if not np.isfinite(x+y):raise ValueError()
  location="wcs"
 except Exception:
  x,y=(nx-1)/2,(ny-1)/2;location="geometry_only"
 yy,xx=np.indices(a.shape);rr=np.hypot(xx-x,yy-y)
 finite=np.isfinite(a);ann=finite&(rr>=6)&(rr<.46*min(nx,ny))
 if ann.sum()<35 or not (finite&(rr<=3)).any():return {"state":"insufficient_pixels"},a
 vals=a[ann];bg=float(np.median(vals))
 noise=float(1.4826*np.median(np.abs(vals-bg)))
 if noise<=0 or not math.isfinite(noise):return {"state":"invalid_noise"},a
 near=finite&(rr<=3);far=finite&(rr>5)
 peak=float(a[near].max());outer=float(a[far].max()) if far.any() else bg
 low=float(vals.min());edge=float(min(x,y,nx-1-x,ny-1-y))
 flags=[]
 if edge<7:flags.append("near_edge")
 if outer>max(3*peak,bg+8*noise):flags.append("bright_neighbour_or_sidelobe")
 if bg-low>8*noise:flags.append("deep_negative_sidelobe_possible")
 return {"state":"audited","WCS":location,"local_noise":round(noise,7),
   "local_snr":round((peak-bg)/noise,2),"center_peak":round(peak,7),
   "outer_peak":round(outer,7),"negative_trough":round(low,7),
   "edge_px":round(edge,2),"flags":";".join(flags)},a
def cutout_one(sid,url,srow,g):
 report={"source":sid,"status":"failed","epochs_expected":len(g),"epochs_audited":0,"url":url}
 results=[];thumb=[]
 try:
  if urlparse(url).hostname!="data.csiro.au":raise ValueError("unexpected remote host")
  b,resolved=get_file(url)
  if not b.startswith(b"SIMPLE"):raise ValueError("remote file not FITS")
  report.update({"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b),"resolved":resolved})
  with fits.open(io.BytesIO(b),memmap=False) as cube:
   images=[h for h in cube[1:] if h.data is not None]
   report["extensions"]=len(images)
   times=[pd.to_datetime(h.header.get("DATE-OBS"),errors="coerce",utc=True) for h in images]
   used=set()
   for _,m in g.sort_values(["time","image_id"]).iterrows():
    q={"source":sid,"image_id":int(m.image_id),"time":str(m.time),
       "forced":flag(m.forced),"catalogue_snr":float(m.snr),"image_state":"unmatched"}
    t=pd.to_datetime(m.time,errors="coerce",utc=True)
    possible=sorted((abs((v-t).total_seconds()),j) for j,v in enumerate(times) if j not in used and pd.notna(t) and pd.notna(v))
    if possible and possible[0][0]<=120:
     dt,j=possible[0];used.add(j)
     h=images[j]
     met,a=pixels(h.data,h.header,float(srow.wavg_ra),float(srow.wavg_dec))
     q.update(met);q.update({"image_state":met["state"],"hdu":j+1,"delta_time_s":round(dt,3)})
     if met["state"]=="audited":
      thumb.append((a,q));report["epochs_audited"]+=1
    results.append(q)
   report["unmatched_extensions"]=len(images)-len(used)
  report["status"]="complete" if report["epochs_audited"]==len(g) else "partial"
  if thumb:
   cols=4;rows=math.ceil(len(thumb)/cols)
   fig,ax=plt.subplots(rows,cols,figsize=(12,rows*2.9),squeeze=False)
   for j,(a,r) in enumerate(thumb):
    p=a[np.isfinite(a)];v0,v1=np.percentile(p,[5,99.5])
    ax.flat[j].imshow(a,origin="lower",cmap="gray",vmin=v0,vmax=max(v1,v0+1e-8))
    ax.flat[j].set_title(f"{r['image_id']} {r['time'][:10]} | SNR {r['local_snr']}",fontsize=8)
    ax.flat[j].set_axis_off()
   for a in ax.flat[len(thumb):]:a.set_axis_off()
   fig.suptitle(f"VAST DR1 {sid} | source FITS-cube native I cutouts",fontsize=10)
   fig.tight_layout()
   (OUT/"montages").mkdir(exist_ok=True)
   fig.savefig(OUT/"montages"/f"{sid}.png",dpi=105)
   plt.close(fig)
 except Exception as e:report["error"]=f"{type(e).__name__}: {str(e)[:170]}"
 return report,results
def gaia_direct(ra,dec,radius):
 """Direct ESAC TAP fallback for CDS outage (the earlier PR 200 route)."""
 endpoint="https://gea.esac.esa.int/tap-server/tap/sync"
 q=("SELECT TOP 25 source_id,ra,dec,parallax,phot_g_mean_mag FROM gaiadr3.gaia_source "
    f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},{radius/3600:.10f}))")
 try:
  r=requests.get(endpoint,params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":q},
     timeout=(12,35),headers={"User-Agent":"Open-Space-LPT-followup/1.1"})
  r.raise_for_status()
  t=pd.read_csv(io.StringIO(r.text),comment="#",dtype={"source_id":str})
  candidates=[]
  for _,x in t.iterrows():
   if val(x.ra) is not None and val(x.dec) is not None:
    sep=angular(ra,dec,float(x.ra),float(x.dec))
    if sep<=radius+.1:
     candidates.append((sep,x))
  if not candidates:return {"status":"no_catalog_match","url":r.url,"catalogue":"Gaia DR3 ESAC TAP"}
  sep,x=min(candidates,key=lambda p:p[0])
  return {"status":"match","arcsec":round(sep,3),"record":{
     "SourceID":str(x.source_id),"RAJ2000":float(x.ra),"DEJ2000":float(x.dec),
     "Plx":val(x.parallax),"Gmag":val(x.phot_g_mean_mag)},
     "url":r.url,"catalogue":"Gaia DR3 ESAC TAP"}
 except Exception as ex:
  return {"status":"unavailable","error":f"{type(ex).__name__}: {str(ex)[:145]}","catalogue":"Gaia DR3 ESAC TAP"}
def parse_tsv(body,ra,dec,radius):
 rows=[];cols=None
 for line in body.splitlines():
  if not line or line.startswith("#"):continue
  q=line.split("\t")
  if cols is None:
   if len(q)>1:cols=[x.strip() for x in q]
   continue
  if len(q)!=len(cols) or all(not x.strip("- ") for x in q):continue
  d=dict(zip(cols,[x.strip() for x in q]))
  for x,y in [("RAJ2000","DEJ2000"),("RAdeg","DEdeg"),("RA_ICRS","DE_ICRS"),("RA","DE")]:
   if val(d.get(x)) is not None and val(d.get(y)) is not None:
    sep=angular(ra,dec,float(d[x]),float(d[y]))
    if sep<=radius+.2:rows.append((sep,d))
    break
 if not rows:return {"status":"no_catalog_match","columns":(cols or [])[:16]}
 sep,d=min(rows,key=lambda x:x[0])
 fields={k:v for k,v in d.items() if v and k in {
  "Name","SourceID","RAJ2000","DEJ2000","RAdeg","DEdeg","Ftot","Fpeak","S_Code",
  "Flag","Plx","e_Plx","Gmag","W1mag","W2mag","gmag","rmag","zmag","objID",
  "PS1","AllWISE","Class","Type","RA_ICRS","DE_ICRS","_r","W1","W2","Stot","S1.4"}}
 return {"status":"match","arcsec":round(sep,3),"record":fields}
def viz(cat,ra,dec,radius):
 """Mirror-aware ASU query; CDS main server outage 9 November 2026."""
 params={"-source":cat,"-c":f"{ra:.8f} {dec:.8f}","-c.rs":str(radius),
         "-out.max":"8","-out.all":"","-out.add":"_r"}
 errors=[]
 for host in ("vizier.cfa.harvard.edu","vizier.idia.ac.za"):
  url=f"https://{host}/viz-bin/asu-tsv"
  try:
   r=requests.get(url,params=params,timeout=(10,19),headers={"User-Agent":"Open-Space-LPT-followup/1.1"})
   r.raise_for_status()
   if r.text.lstrip().startswith("<html") or "Internal Server Error" in r.text[:300]:
    raise ValueError("returned HTML not ASU catalogue rows")
   out=parse_tsv(r.text,ra,dec,radius)
   if out["status"]=="no_catalog_match" and not out["columns"]:
    raise ValueError("no schema returned; not a measured non-detection")
   out.update({"url":r.url,"mirror":host})
   return out
  except Exception as ex:errors.append(f"{host}: {type(ex).__name__}: {str(ex)[:100]}")
 if cat=="I/355/gaiadr3":
  g=gaia_direct(ra,dec,radius)
  g["mirror_errors"]=errors
  return g
 return {"status":"unavailable","catalog":cat,"errors":errors}
def match_one(sid,ra,dec):
 return sid,{name:viz(cat,ra,dec,rad) for name,(cat,rad) in CAT.items()}
def slope(c):
 pts=[]
 for name,nu in FREQ.items():
  item=c.get(name,{})
  if item.get("status")!="match":continue
  f=val(item["record"].get("Ftot"))
  if f is None:f=val(item["record"].get("Fpeak"))
  if f is not None and f>0:pts.append((nu,f))
 if len(pts)<2:return {"state":"insufficient_Stokes_I"}
 pts.sort()
 alpha=math.log(pts[-1][1]/pts[0][1])/math.log(pts[-1][0]/pts[0][0])
 return {"state":"noncontemporaneous_flux_slope","alpha":round(alpha,3),
         "freq_MHz":[pts[0][0],pts[-1][0]],"mJy":[pts[0][1],pts[-1][1]]}
def classify(sid,a,co,img):
 if a and a.get("arcsec",999)<2:return "known_pulsar"
 r=img.get(sid,{})
 if r.get("epochs_audited",0)<r.get("epochs_expected",1):return "unresolved_incomplete_native_images"
 if co.get("gaia",{}).get("status")=="match" and co["gaia"].get("arcsec",999)<2:
  return "stellar_counterpart_flare_unconfirmed"
 if co.get("legacy_dr10_photoz_subset",{}).get("status")=="match":
  return "extragalactic_counterpart_agn_scintillation_possible"
 return "unresolved_not_polarisation_vetted"
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 s,m,ids,on=shortlist()
 atnf={str(x["source"]):x for x in json.loads((BASE/"atnf_receipts.json").read_text())["nearest"]}
 urls,ind=cutout_index(ids);dump("index_receipts.json",ind)
 img={};logs=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for i in range(0,60,8):
   batch=ids[i:i+8]
   futures={pool.submit(cutout_one,x,urls[x],s.loc[x],m[m.source==x]):x for x in batch if x in urls}
   for fut in as_completed(futures):
    sid=futures[fut];img[sid],rr=fut.result();logs+=rr
   for sid in batch:
    if sid not in urls:
     img[sid]={"source":sid,"status":"index_missing","epochs_expected":len(m[m.source==sid]),"epochs_audited":0}
     for _,row in m[m.source==sid].iterrows():
      logs.append({"source":sid,"image_id":int(row.image_id),"time":str(row.time),
                   "forced":flag(row.forced),"catalogue_snr":float(row.snr),"image_state":"index_missing"})
   dump("image_receipts.json",[img[x] for x in ids if x in img])
   print(f"Images {len(img)}/60",flush=True)
 assert len(logs)==693
 pd.DataFrame(logs).sort_values(["source","time","image_id"]).to_csv(OUT/"epoch_image_audit.csv",index=False)
 co={}
 with ThreadPoolExecutor(max_workers=6) as pool:
  for i in range(0,60,10):
   f=[pool.submit(match_one,sid,float(s.loc[sid].wavg_ra),float(s.loc[sid].wavg_dec)) for sid in ids[i:i+10]]
   for done in as_completed(f):
    sid,obj=done.result();co[sid]=obj
   dump("counterpart_receipts.json",co)
   print(f"Surveys {len(co)}/60",flush=True)

 # Independent bounded RACS-low2 circular-polarisation catalogue lookup; never infer upper limits.
 try:
  from racs_v import probe
  vprobe=probe({sid:(float(s.loc[sid].wavg_ra),float(s.loc[sid].wavg_dec)) for sid in ids})
 except Exception as e:
  vprobe={"status":"unavailable","error":f"{type(e).__name__}: {str(e)[:150]}"}
 dump("racs_low2_v_receipts.json",vprobe)
 rows=[]
 for sid in ids:
  c=co.get(sid,{})
  im=[x for x in logs if x["source"]==sid]
  spectrum=slope(c)
  rows.append({"source":sid,"onoff":sid in on,"ra":float(s.loc[sid].wavg_ra),"dec":float(s.loc[sid].wavg_dec),
   "epochs":len(im),"image_audits":sum(x.get("image_state")=="audited" for x in im),
   "image_flags":sum(bool(x.get("flags")) for x in im),
   "psr":atnf[sid]["psr"] if atnf[sid]["arcsec"]<2 else "",
   "class":classify(sid,atnf.get(sid),c,img),
   "racs_alpha":spectrum.get("alpha"),"racs_spectrum_state":spectrum["state"],
   "RACS_V":"not_measured","VAST_V":"not_measured",
   **{k:c.get(k,{}).get("status","unavailable") for k in CAT}})
 pd.DataFrame(rows).to_csv(OUT/"classifications.csv",index=False)
 classes=pd.Series([r["class"] for r in rows]).value_counts().to_dict()
 dump("summary.json",{"status":"conservative_followup_not_discovery","source_count":60,
  "onoff_ids":on,"epoch_measurements":693,
  "native_images":sum(x.get("epochs_audited",0) for x in img.values()),
  "complete_source_images":sum(x.get("status")=="complete" for x in img.values()),
  "survey_queries":60*len(CAT),"classes":classes,
  "stokes_v":"not measured in published VAST source FITS cubes; CASDA AS207/RACS-low2 I,V require separately staged products and leakage controls",
  "periodicity":"not tested: no validated compact polarised counterpart-free target or seconds-scale series",
  "literature":"ATNF v18 prior-identification checked in baseline; individual ATel/ADS check only if a genuine polarised remainder survives",
  "discovery_validated":False,"racs_low2_v_catalogue_probe":vprobe.get("status"),"racs_low2_v_selected_matches":len(vprobe.get("matches",{}))})
 print(json.dumps(classes,indent=2),flush=True)
if __name__=="__main__":
 try:main()
 except Exception as e:
  dump("failure.json",{"error":str(e),"traceback":traceback.format_exc()})
  raise
