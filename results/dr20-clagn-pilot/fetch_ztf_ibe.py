import sys
import urllib.request,urllib.parse,urllib.error,pathlib,hashlib,json
p=pathlib.Path(sys.argv[1]);p.mkdir(parents=True,exist_ok=True)
if any(p.iterdir()): raise SystemExit('Output directory must be empty')
records=[]
for name,lo,hi in [('baseline',59588,59708),('development',60262,60382)]:
 pars={'POS':'213.351140496,53.0908107113','WHERE':f'obsjd >= {lo+2400000.5} AND obsjd <= {hi+2400000.5}','COLUMNS':'obsjd,filtercode,field,ccdid,qid,infobits,filefracday,maglimit,airmass'};u='https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?'+urllib.parse.urlencode(pars);f=p/(name+'.tbl')
 try:
  with urllib.request.urlopen(u,timeout=60) as r:status=r.status;b=r.read()
 except urllib.error.HTTPError as e:status=e.code;b=e.read()
 except Exception as e:records.append(dict(name=name,url=u,error=repr(e)));continue
 f.write_bytes(b);records.append(dict(name=name,url=u,status=status,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),file=f.name));print(name,status,len(b),b[:200])
(p/'requests.json').write_text(json.dumps(records,indent=2)+'\n')
