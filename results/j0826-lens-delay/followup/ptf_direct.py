"""Query PTF epoch detections at J0826+7002 with several indexed services."""
import csv, hashlib, io, json, urllib.parse,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
DIR=Path(__file__).parent/"runner-output"
DIR.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"Open-Space J0826 2026 follow-up"}
def get(url):
 with urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=65) as r:
  return r.read()
def fetch(label,url):
 try:
  raw=get(url)
  ext=".csv" if b"," in raw[:200] else (".xml" if b"<" in raw[:200] else ".tbl")
  (DIR/(label+ext)).write_bytes(raw)
  print("DIRECT",label,url,"BYTES",len(raw),"SHA256",hashlib.sha256(raw).hexdigest(),"HEAD",raw[:350].decode(errors="replace"))
  if b"<" in raw[:200]:
   root=ET.fromstring(raw)
   fields=[x.attrib.get("name") for x in root.findall(".//{*}FIELD")]
   rr=[[t.text for t in x.findall("./{*}TD")] for x in root.findall(".//{*}TR")]
   print("FIELDS",fields[:20],"ROWS",len(rr),"SAMPLE",rr[:2])
  elif b"," in raw[:150]:
   rr=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
   print("ROWS",len(rr),"SAMPLE",json.dumps(rr[:2])[:800])
 except Exception as e: print("DIRECT_ERROR",label,type(e).__name__,str(e)[:600])

if __name__=="__main__":
 base="https://irsa.ipac.caltech.edu"
 pairs=[("ptf_sources",dict(table="ptf_sources",RA=126.5341,DEC=70.0456,SR=.004,format="votable")),
        ("ptf_lightcurves_close",dict(table="ptf_lightcurves",RA=126.5333358,DEC=70.0464092,SR=.00028,format="votable"))]
 for label,args in pairs:
  fetch(label,base+"/SCS?"+urllib.parse.urlencode(args))
 for obj in ("52342070001397","52342070001419"):
  adql="SELECT TOP 250 oid,ra,dec,fid,obsmjd,mag_autocorr,magerr_auto,goodflag,photcalflag,flags,qaflags,limitmag FROM ptf_lightcurves WHERE oid="+obj
  fetch("tap_oid_"+obj,base+"/TAP/sync?"+urllib.parse.urlencode({"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":adql}))
  fetch("gator_oid_"+obj,base+"/cgi-bin/Gator/nph-query?"+urllib.parse.urlencode({"catalog":"ptf_lightcurves","spatial":"None","constraints":"oid="+obj,"outfmt":1}))
