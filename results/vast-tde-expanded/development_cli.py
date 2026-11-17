"""Review-bound VAST development CLI. Never reads heldout rows intentionally."""
import argparse,hashlib,json,platform,subprocess,sys
from pathlib import Path
from importlib.metadata import version
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
SOURCES=(14043503,13894502,14481388)
VERSIONS={'numpy':'2.5.5','scipy':'1.17.1','pyarrow':'24.0.0','astropy':'8.0.2'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def verify_blob_bytes(files,head):
 for name in files:
  blob=subprocess.check_output(['git','-C',str(REPO),'show',f'{head}:{name}'])
  if (REPO/name).read_bytes()!=blob:raise ValueError('worktree/Git blob bytes differ: '+name)

def snapshot():
 if git('status','--porcelain','--','results/vast-tde-expanded'):raise ValueError('packet worktree is not clean')
 head=git('rev-parse','HEAD');files=git('ls-files','results/vast-tde-expanded').splitlines()
 if not files:raise ValueError('packet not tracked')
 verify_blob_bytes(files,head)
 blobs={n:git('rev-parse',f'{head}:{n}') for n in files}
 contents={n:sha(REPO/n) for n in files}
 actual={n:version(n) for n in VERSIONS}
 if actual!=VERSIONS or platform.python_version()!='3.14.0':raise ValueError('dependency version mismatch')
 return {'git_head':head,'git_blobs':blobs,'worktree_sha256':contents,'python':platform.python_version(),'dependencies':actual,'phase':'development','sources':list(SOURCES),'holdout_authorized':False}

def verify_manifest(path,expected,phase):
 if phase!='development':raise ValueError('wrong phase')
 p=Path(path)
 if sha(p)!=expected:raise ValueError('manifest SHA mismatch')
 saved=json.loads(p.read_text(encoding='utf-8'));now=snapshot()
 if saved!=now:raise ValueError('code, packet, HEAD, or environment mismatch')
 return saved

def d_path(p):
 p=Path(p).resolve()
 if p.drive.upper()!='D:':raise ValueError('output must be on D:')
 return p

def run(args):
 if not args.authorize_development:raise ValueError('explicit development authorization required')
 if args.phase!='development':raise ValueError('wrong phase')
 manifest=verify_manifest(args.manifest,args.manifest_sha256,args.phase)
 output=d_path(args.output)
 if output.exists():raise FileExistsError('output directory must be new')
 # Hash external inputs before any code can decode radio flux.
 from read_development import digest
 if digest(Path(args.images))!='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195':raise ValueError('image metadata digest mismatch')
 from packet_guard import guard
 packet=guard()
 if digest(Path(args.parquet))!=packet['measurement_parquet_sha256']:raise ValueError('parquet digest mismatch')
 output.mkdir(parents=True)
 from read_development import read_development
 from conditional_pipeline import develop,freeze_development
 import csv
 import numpy as np
 import traceback
 controls=list(csv.DictReader((ROOT/'control-metadata.csv').open(newline='',encoding='utf-8')))
 report={'manifest_sha256':args.manifest_sha256,'git_head':manifest['git_head'],'sources':{},'conditional_only':True}
 for source in SOURCES:
  dropped=[]
  try:
   rows,dropped=read_development(source,args.parquet,args.images)
   ids={int(r['control_source']) for r in controls if int(r['target_source'])==source}
   result=develop(source,packet,rows,ids);result['drops']=dropped+result['drops']
   path=output/f'{source}-development.json'
   digest_value=freeze_development(path,result,args.manifest_sha256,manifest['worktree_sha256']['results/vast-tde-expanded/conditional_pipeline.py'])
   report['sources'][str(source)]={'status':'qualified-for-review','file':path.name,'sha256':digest_value,'input_rows':len(rows)+len(dropped),'drops':dropped}
  except Exception as exc:
   report['sources'][str(source)]={'status':'inconclusive','error_type':type(exc).__name__,'error':str(exc),'drops':dropped}
 with (output/'summary.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
 return report

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--make-manifest',type=Path);ap.add_argument('--manifest',type=Path);ap.add_argument('--manifest-sha256');ap.add_argument('--phase',default='development');ap.add_argument('--authorize-development',action='store_true');ap.add_argument('--parquet',type=Path);ap.add_argument('--images',type=Path);ap.add_argument('--output',type=Path);a=ap.parse_args()
 if a.make_manifest:
  target=d_path(a.make_manifest);value=snapshot();raw=(json.dumps(value,indent=2,sort_keys=True)+'\n').encode()
  with target.open('xb') as f:f.write(raw)
  print('manifest SHA-256',hashlib.sha256(raw).hexdigest())
 else:
  if not all((a.manifest,a.manifest_sha256,a.parquet,a.images,a.output)):ap.error('manifest, manifest-sha256, parquet, images and output required')
  print(json.dumps(run(a),indent=2))
