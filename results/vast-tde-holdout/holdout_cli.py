"""Sealed holdout execution. No data access without separate reviewed authorization.
Parquet predicates limit returned rows, not physical column-chunk decoding.
"""
import argparse,csv,hashlib,json,platform,subprocess,sys
from importlib.metadata import version
from pathlib import Path
ROOT=Path(__file__).resolve().parent;REPO=ROOT.parents[1];BASE=ROOT.parent/'vast-tde-expanded';sys.path.insert(0,str(BASE))
OLD_HEAD='f323ab5e2b4c0b9d0227c581884c353e6948e39a'
OLD_MANIFEST_SHA='21ccdd9daeff9cb52c03eb0c2a469b005b1a63eb4c65c9a8e6a727d4a0cd8a5d'
SUMMARY_SHA='a851263393a0f59609fa14da6443e9c9bb6b751b94c4821b96da9af4f6e4bb60'
DEV_SHA={'14043503':'40d07ba0dc985afd396e166c6529fc92a9ede3f5b4dbca5c028fff22a501901f','13894502':'efd15b30c99b7adacddf9868348a5312530733240da7dfa587492f21df336941','14481388':'b2ffea8dce960b8aa1cbd86b88f8d160f2bb373562f4e9e2f6ecd13cca2d036e'}
VERSIONS={'numpy':'2.5.5','scipy':'1.17.1','pyarrow':'24.0.0','astropy':'8.0.2'}
IMAGE_SHA='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195'
PARQUET_SHA='a3af873784b4b3d5a727d6d50890a48a5fd1bbc41ebe035106ce83d47cfb9ca6'

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def snapshot():
 if git('status','--porcelain','--','results/vast-tde-holdout','results/vast-tde-expanded'):raise ValueError('dirty packet')
 head=git('rev-parse','HEAD');files=git('ls-files','results/vast-tde-holdout','results/vast-tde-expanded').splitlines();blobs={};work={}
 for name in files:
  b=subprocess.check_output(['git','-C',str(REPO),'show',f'{head}:{name}'])
  if (REPO/name).read_bytes()!=b:raise ValueError('worktree/blob mismatch '+name)
  blobs[name]=git('rev-parse',f'{head}:{name}');work[name]=hashlib.sha256(b).hexdigest()
 actual={n:version(n) for n in VERSIONS}
 if actual!=VERSIONS or platform.python_version()!='3.14.0':raise ValueError('dependency mismatch')
 return dict(git_head=head,old_head=OLD_HEAD,old_manifest_sha256=OLD_MANIFEST_SHA,old_summary_sha256=SUMMARY_SHA,dev_sha256=DEV_SHA,git_blobs=blobs,worktree_sha256=work,python=platform.python_version(),dependencies=actual,phase='holdout',authorization_required=True)

def verify_old(old_manifest,dev_dir):
 if sha(old_manifest)!=OLD_MANIFEST_SHA:raise ValueError('original manifest mismatch')
 old=json.loads(Path(old_manifest).read_text(encoding='utf-8'))
 if old['git_head']!=OLD_HEAD:raise ValueError('original head mismatch')
 for name,digest in old['worktree_sha256'].items():
  if sha(REPO/name)!=digest:raise ValueError('original code or packet changed: '+name)
 summary=Path(dev_dir)/'summary.json'
 if sha(summary)!=SUMMARY_SHA:raise ValueError('development summary mismatch')
 d=json.loads(summary.read_text(encoding='utf-8'))
 if d['manifest_sha256']!=OLD_MANIFEST_SHA or d['git_head']!=OLD_HEAD or set(d['sources'])!=set(DEV_SHA):raise ValueError('development provenance mismatch')
 result={}
 for sid,expected in DEV_SHA.items():
  entry=d['sources'][sid]
  if entry['status']!='qualified-for-review' or entry['sha256']!=expected:raise ValueError('development summary entry mismatch')
  f=Path(dev_dir)/entry['file']
  if sha(f)!=expected:raise ValueError('development file mismatch '+sid)
  v=json.loads(f.read_text(encoding='utf-8'))
  if v['packet_sha256']!=OLD_MANIFEST_SHA or v['code_sha256']!=old['worktree_sha256']['results/vast-tde-expanded/conditional_pipeline.py']:raise ValueError('development code provenance mismatch')
  result[sid]=v['result']
 return result

def verify_review(manifest,manifest_sha,authorization,auth_sha,phase):
 if phase!='holdout':raise ValueError('wrong phase')
 if sha(manifest)!=manifest_sha:raise ValueError('holdout manifest mismatch')
 m=json.loads(Path(manifest).read_text(encoding='utf-8'))
 if m!=snapshot():raise ValueError('holdout code or packet mismatch')
 if sha(authorization)!=auth_sha:raise ValueError('authorization digest mismatch')
 a=json.loads(Path(authorization).read_text(encoding='utf-8'))
 if a.get('phase')!='holdout' or a.get('authorized') is not True or a.get('git_head')!=m['git_head'] or a.get('manifest_sha256')!=manifest_sha or a.get('development_summary_sha256')!=SUMMARY_SHA:raise ValueError('authorization does not match reviewed packet')
 return m

def holdout_rows(source,frozen,packet,parquet,images):
 # Hash all files before reading any Parquet science columns.
 if sha(images)!=IMAGE_SHA or sha(parquet)!=PARQUET_SHA:raise ValueError('input digest mismatch')
 import pyarrow.parquet as pq
 from conditional_pipeline import phase_rows
 item=packet['sources'][str(source)];groups=item['holdout']
 if groups!=frozen['holdout_images'] or len(groups)!=2:raise ValueError('holdout group mismatch')
 ids={i for g in groups for i in g};dev={i for g in item['development'] for i in g}
 if ids&dev:raise ValueError('phase overlap')
 controls={int(k) for k in frozen['control_baselines']};sources=controls|{source}
 columns=['source','image_id','flux_peak','flux_peak_err','local_rms','forced','flag_c4','has_siblings']
 tab=pq.read_table(parquet,columns=columns,filters=[('source','in',list(sources)),('image_id','in',list(ids))])
 if not set(tab.column('source').to_pylist())<=sources or not set(tab.column('image_id').to_pylist())<=ids:raise ValueError('forbidden returned rows')
 meta={}
 with Path(images).open(newline='',encoding='utf-8') as f:
  for r in csv.DictReader(f):
   if int(r['image_id']) in ids:meta[int(r['image_id'])]=(r['beam_bmaj'],r['beam_bmin'])
 if set(meta)!=ids:raise ValueError('missing image metadata')
 rows=tab.to_pylist()
 for r in rows:r['beam_bmaj'],r['beam_bmin']=meta[int(r['image_id'])]
 return phase_rows(rows,ids,sources)

def run(args):
 if not args.authorize_holdout:raise ValueError('separate holdout authorization required')
 m=verify_review(args.manifest,args.manifest_sha256,args.authorization,args.authorization_sha256,args.phase)
 baseline=verify_old(args.old_manifest,args.development_dir)
 from packet_guard import guard
 packet=guard();out=Path(args.output).resolve()
 if out.drive.upper()!='D:' or out.exists():raise ValueError('new D: output directory required')
 if sha(args.images)!=IMAGE_SHA or sha(args.parquet)!=PARQUET_SHA:raise ValueError('input digest mismatch')
 out.mkdir(parents=True)
 from conditional_pipeline import evaluate_holdout
 report={'git_head':m['git_head'],'manifest_sha256':args.manifest_sha256,'authorization_sha256':args.authorization_sha256,'old_summary_sha256':SUMMARY_SHA,'scope':'conditional standard-release repeatability; not physical significance','sources':{}}
 for sid in DEV_SHA:
  drops=[]
  try:
   rows,drops=holdout_rows(int(sid),baseline[sid],packet,args.parquet,args.images)
   value=evaluate_holdout(baseline[sid],rows);value['drops']=drops+value['drops'];report['sources'][sid]={'status':'measured',**value}
  except Exception as exc:report['sources'][sid]={'status':'inconclusive','error_type':type(exc).__name__,'error':str(exc),'drops':drops}
 with (out/'holdout-result.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
 return report

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--make-manifest',type=Path);ap.add_argument('--manifest');ap.add_argument('--manifest-sha256');ap.add_argument('--authorization');ap.add_argument('--authorization-sha256');ap.add_argument('--phase',default='holdout');ap.add_argument('--authorize-holdout',action='store_true');ap.add_argument('--old-manifest');ap.add_argument('--development-dir');ap.add_argument('--parquet');ap.add_argument('--images');ap.add_argument('--output');a=ap.parse_args()
 if a.make_manifest:
  p=a.make_manifest.resolve()
  if p.drive.upper()!='D:':raise ValueError('manifest must be on D:')
  raw=(json.dumps(snapshot(),sort_keys=True,indent=2)+'\n').encode()
  with p.open('xb') as f:f.write(raw)
  print('manifest SHA-256',hashlib.sha256(raw).hexdigest())
 else:
  if not all((a.manifest,a.manifest_sha256,a.authorization,a.authorization_sha256,a.old_manifest,a.development_dir,a.parquet,a.images,a.output)):ap.error('all holdout inputs required')
  print(json.dumps(run(a),indent=2))
