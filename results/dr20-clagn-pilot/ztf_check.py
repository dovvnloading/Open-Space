"""Check exact IRSA response bytes and frozen positional/quality coverage."""
import csv,gzip,hashlib,io,json,math
from pathlib import Path
p=Path(__file__).resolve().parent/'ztf';rec=json.loads((p/'requests.json').read_text());results=[]
for x in rec:
 b=gzip.decompress((p/(x['file']+'.gz')).read_bytes());assert hashlib.sha256(b).hexdigest()==x['sha256'];rows=list(csv.DictReader(io.StringIO(b.decode())));target=[];controls=[]
 for r in rows:
  d=3600*math.hypot((float(r['ra'])-213.351140496)*math.cos(math.radians(53.0908107113)),float(r['dec'])-53.0908107113)
  if int(r['catflags'])!=0 or not math.isfinite(float(r['magerr'])) or float(r['magerr'])<=0:continue
  if d<=2:target.append(r)
  if 20<=d<=60:controls.append(r)
 results.append(dict(file=x['file'],raw_rows=len(rows),good_target_rows=len(target),good_control_rows=len(controls),good_control_ids=sorted(set(r['oid'] for r in controls))))
print(json.dumps(results,indent=2));assert all(x['good_target_rows']==0 for x in results)
