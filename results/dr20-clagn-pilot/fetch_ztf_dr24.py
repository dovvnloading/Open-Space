"""Bounded explicit DR24 query control; never overwrites existing receipts."""
import urllib.request,urllib.parse,urllib.error,pathlib,hashlib,json,sys,concurrent.futures
p=pathlib.Path(sys.argv[1]);p.mkdir(parents=True,exist_ok=True)
if any(p.iterdir()):raise SystemExit('output directory must be empty')
jobs=[]
for ra,label in [('213.351140496','normal'),('-146.648859504','wrapped')]:
 for band in ('g','r'):
  for name,t in [('all',None),('baseline','59588 59708'),('development','60262 60382')]:
   if label=='wrapped' and t is not None:continue
   pars={'POS':f'CIRCLE {ra} 53.0908107113 0.001388889','BANDNAME':band,'COLLECTION':'ztf_dr24','FORMAT':'CSV'}
   if t:pars['TIME']=t
   jobs.append((f'{label}-{band}-{name}',pars))
def get(job):
 name,pars=job;u='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?'+urllib.parse.urlencode(pars);path=p/(name+'.response')
 try:
  with urllib.request.urlopen(u,timeout=70) as r:status=r.status;headers=dict(r.headers);body=r.read()
 except urllib.error.HTTPError as e:status=e.code;headers=dict(e.headers);body=e.read()
 except Exception as e:return dict(name=name,url=u,error=repr(e))
 path.write_bytes(body);return dict(name=name,url=u,status=status,content_type=headers.get('Content-Type'),bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),file=path.name)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:out=list(ex.map(get,jobs))
(p/'requests.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
