import sys
import urllib.request,urllib.parse,pathlib,hashlib,json
p=pathlib.Path(sys.argv[1]);p.mkdir(parents=True,exist_ok=True)
if any(p.iterdir()): raise SystemExit('Output directory must be empty')
rec=[]
for band in ['g','r']:
 for label,lo,hi in [('baseline',59588,59708),('development',60262,60382)]:
  args=dict(POS='CIRCLE 213.351140496 53.0908107113 0.0166667',BANDNAME=band,TIME=f'{lo} {hi}',FORMAT='CSV');u='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?'+urllib.parse.urlencode(args);path=p/f'{band}-{label}.csv'
  try:
   with urllib.request.urlopen(u,timeout=90) as r:b=r.read();status=r.status;headers=dict(r.headers)
   path.write_bytes(b);rec.append(dict(file=path.name,url=u,status=status,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),content_type=headers.get('Content-Type')));print(band,label,len(b),b[:250],flush=True)
  except Exception as e:rec.append(dict(file=path.name,url=u,error=repr(e)));print(band,label,repr(e),flush=True)
(p/'requests.json').write_text(json.dumps(rec,indent=2)+'\n')
