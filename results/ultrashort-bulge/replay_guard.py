"""Fail-closed primitives for immutable microlensing replay campaigns."""
import csv,hashlib,json,math,os,platform,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
FROZEN={'kmt/pySIS-2024.tsv':'863f215cf1a908c19148c8a4de15cd94d13474fe8b55fa76712ee23517e7391f','kmt/pySIS-2025.tsv':'d6e3edfdf28b29c43b62df29a4e4b240db9c2e1cf4046249ee1fddde1b2f23ed','ogle/catalog.csv':'7e9be993664be633f5225d72726ffe2bcafc557a4b7920c81e88f09b9d2edd1c','kmt/catalog-2024.tsv':'77187fe914f5da33a0ae8169ebab7dcc3fb3f3977f2231e302d9b1fedc6a9bb5','kmt/catalog-2025.tsv':'e2947ff7578590cd48192fdbb5fdc236bf17d9207835f50b7af1b5ba45bbb791','tables/cadence.csv.gz':'9d2d56ceebb1512da151edb1b23acbae02d7f3a02decafb68e24cd0ec31b2bdb','tables/peak-comparisons.json.gz':'3e87ce4e395e73a70b39b0f2d2b5e8247e5fa3c0d046b7fc1a6c04c426376567','other/gaia.csv':'52a394f0df36e6baf6b3c31909dd43d9c5abf0a2a994fa5c28b420580fdde80a','other/gaia_index_receipt.json':'3b86fabcf9470a26ea89275641b5dce8b910730c538d15b58861171e8330a8f9','ogle/tables.json':'ee71ddb912e1d60b8f2eee228f58926e238ee93c78f83904fb726b4a3da0811a','ogle/photometry.csv':'4cdf77ca3ca3cc98bf180084d8097b5838f95a6e10dc85ea3ce2f8df49c6dd6a','kmt/index-receipts.json':'6c60b82229d10a228bc0a2e737312f00a1f1dea1ff0cef1cb13eff5c57afd37b','other/index_receipts.csv':'2a632070d11636b52a06717a8ca47b40119b122ccd40e4e7bcd692395b8ae075','other/photometry.csv':'ad458410df492307f0b7beb73645142118ab07f5403c5e69d05dff749dfeb6d0','tables/grid-fits.csv.gz':'dafc904395f8a0771c2beda9c820b1f3cef216a95b734035a0b841b6ece234af'}
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def external(path):
 p=Path(path)
 if not p.is_absolute():raise ValueError('absolute D: path required')
 p=p.resolve()
 allowed=(p.drive.upper()=='D:') if os.name=='nt' else (p==Path('/mnt/d') or Path('/mnt/d') in p.parents)
 if not allowed or p==REPO or REPO in p.parents:raise ValueError('external D: path required')
 return p
def fresh(path):
 p=external(path)
 if p.exists():raise FileExistsError('new exclusive output required: '+str(p))
 p.mkdir(parents=True,exist_ok=False);return p
def verified(path,expected):
 p=Path(path)
 if not p.is_file() or sha(p)!=expected:raise ValueError('frozen input mismatch: '+str(p))
 return p
def identity(path,event,cap,code_hashes,runtime,parameters=None):
 return {'event':event,'input_sha256':sha(path),'cap':bool(cap),'code_sha256':dict(code_hashes),'runtime':runtime,'parameters':parameters}
def check_checkpoint(record,expected):
 if record.get('_identity')!=expected:raise ValueError('stale or unbound checkpoint')
 return record

def verify_receipt(name):
 return verified(ROOT/name,FROZEN[name])

def fresh_for_inputs(out,*inputs):
 dest=external(out)
 for item in inputs:
  src=external(item)
  if dest==src or dest in src.parents or src in dest.parents:raise ValueError('input/output overlap')
 return fresh(dest)

def frozen_parameters(cap):
 for name in ('ogle/catalog.csv','kmt/catalog-2024.tsv','kmt/catalog-2025.tsv'):verify_receipt(name)
 def number(value):
  try:
   x=float(value);return x if math.isfinite(x) else None
  except (TypeError,ValueError):return None
 result={}
 if not cap:
  with (ROOT/'ogle/catalog.csv').open(newline='',encoding='utf-8-sig') as f:
   for r in csv.DictReader(f):result['OGLE-'+r['event']]={'t0':number(r['t0_hjd']),'te':number(r['tau_days'])}
 for year in (2024,2025):
  with (ROOT/'kmt'/f'catalog-{year}.tsv').open(newline='',encoding='utf-8-sig') as f:
   for r in csv.DictReader(f,delimiter='\t'):
    t=number(r['t_0']);te=number(r['t_E'])
    result[r['Event']]={'t0':2400000+t if t is not None and t>0 else None,'te':te if te is not None and te>0 else None}
 return result

def parameters_digest(parameters):
 return hashlib.sha256(json.dumps(parameters,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def validate_job_membership(jobs,expected):
 names=[j[0] for j,_ in jobs]
 if len(names)!=len(set(names)):raise ValueError('duplicate event IDs in jobs')
 if set(names)!=set(expected):raise ValueError('incomplete or changed frozen event membership')

def engine_assets(vbm):
 package=Path(vbm.__file__).resolve()
 module=sys.modules[vbm.VBMicrolensing.__module__]
 native=Path(module.__file__).resolve()
 table=package.parent/'data'/'ESPL.tbl'
 if not all(p.is_file() for p in (package,native,table)):raise ValueError('missing microlensing engine asset')
 return (package,native,table)

CAMPAIGN_SCHEMA='ultrashort-grid-campaign-v1'
RUNTIME_KEYS=('python','platform','numpy','vbm')
ENGINE_ROLES=('vbm_package','vbm_native','espl_table','grid_library')

def source_assets():
 return {'run_grid':ROOT/'run_grid.py','fit':ROOT/'fit.py','grid_source':ROOT/'grid.cpp','kmt_flux':ROOT/'kmt/flux.py','guard':ROOT/'replay_guard.py','ogle_catalog':ROOT/'ogle/catalog.csv','kmt_catalog_2024':ROOT/'kmt/catalog-2024.tsv','kmt_catalog_2025':ROOT/'kmt/catalog-2025.tsv','lens':ROOT.parent/'ogle-dark-lenses/lens.py'}

def validate_campaign(meta):
 if meta.get('schema')!=CAMPAIGN_SCHEMA or meta.get('model')!='frozen grid' or not isinstance(meta.get('cap'),bool):raise ValueError('unsupported campaign schema/model')
 runtime=meta.get('runtime')
 if not isinstance(runtime,dict) or set(runtime)!=set(RUNTIME_KEYS) or any(not isinstance(runtime[k],str) or not runtime[k].strip() or runtime[k].strip().lower() in ('unknown','none') for k in RUNTIME_KEYS):raise ValueError('missing campaign runtime identity')
 roles=meta.get('asset_roles');hashes=meta.get('code_sha256');expected=source_assets()
 if not isinstance(roles,dict) or set(roles)!=set(expected)|set(ENGINE_ROLES):raise ValueError('missing campaign asset role')
 if not isinstance(hashes,dict) or not hashes or len(set(roles.values()))!=len(roles) or set(hashes)!=set(roles.values()):raise ValueError('invalid campaign asset mapping')
 if any(not isinstance(v,str) or re.fullmatch('[0-9a-f]{64}',v) is None for v in hashes.values()):raise ValueError('invalid campaign asset digest')
 for role,path in expected.items():
  if Path(roles[role]).resolve()!=path.resolve():raise ValueError('campaign source identity changed: '+role)
 package=Path(roles['vbm_package']).resolve();native=Path(roles['vbm_native']).resolve();table=Path(roles['espl_table']).resolve()
 if package.name!='__init__.py' or package.parent.name!='VBMicrolensing' or native.parent!=package.parent or not native.name.startswith('VBMicrolensing') or native.suffix.lower() not in ('.so','.pyd','.dll') or table!=package.parent/'data'/'ESPL.tbl':raise ValueError('invalid numerical engine identity')
 if Path(roles['grid_library']).suffix.lower() not in ('.so','.dll','.dylib'):raise ValueError('invalid native grid library')
 for name,digest in hashes.items():
  if sha(name)!=digest:raise ValueError('checkpoint asset changed: '+name)
 return meta
