#!/usr/bin/env python3
"""Direct public optical/IR/X-ray archive cone searches for VAST DR1; fail closed."""
import io, math
import pandas as pd
import requests
from astropy.table import Table
from astropy.io.votable import parse_single_table
AGENT={"User-Agent":"Open-Space-native-VAST-DR1-followup/2.0"}
def arcsec(ra,dec,a,b):
 dx=math.radians((a-ra+180)%360-180)
 d1,d2=math.radians(dec),math.radians(b)
 z=math.sin((d2-d1)/2)**2+math.cos(d1)*math.cos(d2)*math.sin(dx/2)**2
 return 7200*math.degrees(math.asin(min(1,math.sqrt(max(0,z)))))
def number(x):
 try:
  y=float(x)
  return y if math.isfinite(y) else None
 except (TypeError,ValueError):return None
def pick(table,ra,dec,radius,poskeys):
 """Require coordinate columns, not a nearby mismatched source or an error page."""
 if not set(poskeys)<=set(table.colnames):
  raise ValueError("no expected coordinate fields in actual catalogue")
 if len(table)>10000:raise ValueError("cone returned unexpectedly many rows")
 match=[]
 for row in table:
  a,b=number(row[poskeys[0]]),number(row[poskeys[1]])
  if a is None or b is None:continue
  delta=arcsec(ra,dec,a,b)
  if delta<=radius+.1:match.append((delta,row))
 if not match:return {"status":"no_catalog_match","rows_in_cone":len(table)}
 delta,row=min(match,key=lambda z:z[0])
 allowed=("ra","dec","raMean","decMean","objID","nDetections","gMeanPSFMag","rMeanPSFMag",
   "type","flux_g","flux_r","flux_z","w1flux","w2flux","w1mpro","w2mpro",
   "unwise_objid","designation","source_id","srcid","name","_2sxps","id")
 record={}
 for k in table.colnames:
  if k in allowed:
   v=row[k]
   if number(v) is not None:record[k]=float(v)
   else:record[k]=str(v)[:120]
 return {"status":"match","arcsec":round(delta,3),"rows_in_cone":len(table),"record":record}
def csv_table(r):
 if not r.text or "<html" in r.text[:1000].lower():raise ValueError("not machine-readable CSV")
 df=pd.read_csv(io.StringIO(r.text),comment="#")
 return Table.from_pandas(df)
def vo_table(r):
 if not r.content.lstrip().startswith((b"<?xml",b"<VOTABLE",b"<votable")):
  raise ValueError("not machine-readable VOTable")
 t=parse_single_table(io.BytesIO(r.content))
 table=t.to_table(use_names_over_ids=True)
 # Reject a server-side QUERY_STATUS=ERROR masquerading as empty matches
 return table
def attempt(name,endpoint,params,ra,dec,radius,kind,pos):
 out={"catalogue":name,"status":"unavailable","endpoint":endpoint,"search_arcsec":radius}
 try:
  r=requests.get(endpoint,params=params,timeout=(10,28),headers=AGENT)
  out.update({"http":r.status_code,"url":r.url})
  r.raise_for_status()
  if len(r.content)>6_000_000:raise ValueError("oversized cone response")
  tab=vo_table(r) if kind=="vot" else csv_table(r)
  out.update(pick(tab,ra,dec,radius,pos))
  out["schema"]=[str(x) for x in tab.colnames[:20]]
 except Exception as ex:
  out.update({"status":"unavailable","error":f"{type(ex).__name__}: {str(ex)[:170]}"})
 return out
def legacy(ra,dec):
 # NOIRLab Astro Data Lab: complete Legacy Surveys DR10 tractor catalogue
 rad=3
 q=(f"SELECT TOP 16 ra,dec,type,flux_g,flux_r,flux_z FROM ls_dr10.tractor "
    f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},{rad/3600:.10f}))")
 url="https://datalab.noirlab.edu/tap/sync"
 return attempt("NOIRLab Legacy Surveys DR10 tractor",url,
   {"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":q},
    ra,dec,rad,"csv",("ra","dec"))
def unwise(ra,dec):
 rad=3
 return attempt("NASA/IPAC unWISE 2019", "https://irsa.ipac.caltech.edu/SCS",
   {"table":"unwise_2019","RA":f"{ra:.8f}","DEC":f"{dec:.8f}","SR":rad/3600},
   ra,dec,rad,"vot",("ra","dec"))
def panstarrs(ra,dec):
 rad=3
 return attempt("MAST Pan-STARRS DR1 mean", "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr1/mean.csv",
   {"ra":f"{ra:.8f}","dec":f"{dec:.8f}","radius":rad/3600,"pagesize":15,
    "columns":"[objID,raMean,decMean,nDetections,gMeanPSFMag,rMeanPSFMag]"},
   ra,dec,rad,"csv",("raMean","decMean"))
def swift(ra,dec):
 rad=12
 return attempt("NASA HEASARC Swift 2SXPS",
   "https://heasarc.gsfc.nasa.gov/xamin/vo/cone",
   {"table":"swift2sxps","RA":f"{ra:.8f}","DEC":f"{dec:.8f}","SR":rad/3600},
   ra,dec,rad,"vot",("ra","dec"))
QUERIES={"legacy_dr10":legacy,"unwise":unwise,"panstarrs":panstarrs,"swift":swift}
def run(name,ra,dec):
 return QUERIES[name](ra,dec)
