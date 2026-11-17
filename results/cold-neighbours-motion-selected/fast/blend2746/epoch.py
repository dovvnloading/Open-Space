from pathlib import Path
import csv,math,os,hashlib,json
from urllib.request import urlopen
import numpy as np
from astropy.io import fits
R=Path(__file__).resolve().parent
r=next(x for x in csv.DictReader((R.parent/"tracks.csv").open()) if x['source_name']=='track2746p666_1')
ra,dec,mua,mud=map(float,(r['ra2015_deg'],r['dec2015_deg'],r['pmra'],r['pmdec']))
ids=set(r['ids'].split())
out=[]
m=json.loads((R.parent/"manifest.json").read_text())
met=next(field["urls"] for field in m["fields"] if field["tile"]=="2746p666")
cache=Path(os.environ.get("UNTIMELY_CACHE",R/"cache"));cache.mkdir(exist_ok=True,parents=True)
for entry in met:
 p=cache/entry["url"].rsplit("/",1)[-1]
 if not p.exists():p.write_bytes(urlopen(entry["url"],timeout=150).read())
 if hashlib.sha256(p.read_bytes()).hexdigest()!=entry["sha256"]:raise ValueError(f"hash mismatch: {p.name}")
 tab=fits.getdata(p,1);t=(tab[0]['MJDMEAN']-57170)/365.25
 x=(tab['ra']-ra)*3600*math.cos(math.radians(dec));y=(tab['dec']-dec)*3600
 for i in np.flatnonzero(np.hypot(x,y)<25):
  a=tab[i]
  out.append(dict(id=a['unwise_detid'],used=int(a['unwise_detid'] in ids),mjd=a['MJDMEAN'],fwd=a['FORWARD'],x=x[i],y=y[i],mag=22.5-2.5*np.log10(a['flux']) if a['flux']>0 else np.nan,snr=a['flux']/a['dflux'],dx=a['dx'],dy=a['dy'],nm=a['nm'],fracflux=a['fracflux'],flags=a['flags_unwise'],sep=np.hypot(x[i]-mua*t,y[i]-mud*t)))
with (R/'neighbors.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
for a in sorted((v for v in out if v['used']),key=lambda x:x['mjd']):
 print(round(a['mjd']),a['fwd'],round(a['x'],2),round(a['y'],2),round(a['mag'],2),round(a['snr'],1),round(a['fracflux'],2))
