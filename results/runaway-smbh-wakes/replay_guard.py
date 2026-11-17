"""Fail-closed filesystem and frozen cutout checks for research replays."""
import csv,hashlib,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]

def absolute_external(path):
 p=Path(path)
 if not p.is_absolute():raise ValueError('absolute D: path required')
 p=p.resolve()
 if p.drive.upper()!='D:':raise ValueError('absolute D: path required')
 if p==REPO or REPO in p.parents:raise ValueError('output must be outside checkout')
 return p

def fresh_dirs(*paths):
 paths=[absolute_external(p) for p in paths]
 if any(p.exists() for p in paths):raise FileExistsError('output/cache exists; use a new campaign')
 keys=[os.path.normcase(str(p)) for p in paths]
 if len(set(keys))!=len(keys):raise ValueError('output paths collide')
 if any(a in b.parents or b in a.parents for i,a in enumerate(paths) for b in paths[i+1:]):raise ValueError('output paths overlap')
 for p in paths:p.mkdir(parents=True,exist_ok=False)
 return paths

def fresh_output(input_dir,out):
 inp=absolute_external(input_dir);dest=absolute_external(out)
 if inp==dest or inp in dest.parents or dest in inp.parents:raise ValueError('input/output overlap')
 return fresh_dirs(dest)[0]

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def frozen_tiles(tiles_path):
 if sha(tiles_path)!='bfd3724847646662e651363d4624728b69cb9a50577548320e5d7fddf89494b6':raise ValueError('tile list differs from frozen selection')
 with Path(tiles_path).open(newline='',encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 return [r['tile'] for r in rows if r['field']!='LDN1641']

def validate_cutouts(input_dir,tiles,receipt_path=None):
 input_dir=absolute_external(input_dir)
 if receipt_path is None:
  receipt_path=ROOT/'euclid/q1_cutouts.csv'
  if sha(receipt_path)!='a987e71517e3513dabfa6b6b3dae46ccaf1b69cbb3ab6356bb0dd98c557c1e7b':raise ValueError('frozen receipt changed')
 with Path(receipt_path).open(newline='',encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 for tile in tiles:
  for kind in ('VIS','rms'):
   found=[r for r in rows if r['tile']==tile and r['product']==kind]
   if len(found)!=1:raise ValueError('missing or duplicate receipt '+tile+' '+kind)
   p=input_dir/(tile+'_'+kind+'.fits');r=found[0]
   if not p.is_file() or p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:raise ValueError('cutout bytes mismatch '+tile+' '+kind)
 return True

def seal(directory,kind):
 directory=absolute_external(directory)
 files={p.name:sha(p) for p in sorted(directory.iterdir()) if p.is_file()}
 code={p.relative_to(ROOT).as_posix():sha(p) for p in [ROOT/'replay_guard.py',*sorted((ROOT/'detector').glob('wake_*.py')),ROOT/'detector/screen_centers.py',ROOT/'detector/calibrate_centers.py']}
 record={'kind':kind,'files':files,'code_sha256':code,'tiles_sha256':'bfd3724847646662e651363d4624728b69cb9a50577548320e5d7fddf89494b6','receipts_sha256':'a987e71517e3513dabfa6b6b3dae46ccaf1b69cbb3ab6356bb0dd98c557c1e7b'}
 with (directory/'checkpoint_manifest.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(record,f,indent=2);f.write('\n')

def verify_seal(directory,kind):
 directory=absolute_external(directory);p=directory/'checkpoint_manifest.json';record=json.loads(p.read_bytes())
 if record['kind']!=kind or record['tiles_sha256']!='bfd3724847646662e651363d4624728b69cb9a50577548320e5d7fddf89494b6' or record['receipts_sha256']!='a987e71517e3513dabfa6b6b3dae46ccaf1b69cbb3ab6356bb0dd98c557c1e7b':raise ValueError('wrong checkpoint identity')
 actual={p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!='checkpoint_manifest.json'}
 if not actual or actual!=record['files']:raise ValueError('checkpoint bytes or membership changed')
 for name,expected in record['code_sha256'].items():
  if sha(ROOT/name)!=expected:raise ValueError('checkpoint code changed: '+name)
 return actual
