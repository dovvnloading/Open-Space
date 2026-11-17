"""Post-unblinding compact projections of only the authorized source/image rows."""
import argparse,csv,gzip,hashlib,json,sys
from pathlib import Path
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'vast-tde-expanded'))
from packet_guard import guard
from conditional_pipeline import phase_rows
COLUMNS=['source','image_id','flux_peak','flux_peak_err','local_rms','forced','flag_c4','has_siblings']
DATA_SHA='a3af873784b4b3d5a727d6d50890a48a5fd1bbc41ebe035106ce83d47cfb9ca6'
IMAGE_SHA='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195'
DEV_SHA={'14043503':'40d07ba0dc985afd396e166c6529fc92a9ede3f5b4dbca5c028fff22a501901f','13894502':'efd15b30c99b7adacddf9868348a5312530733240da7dfa587492f21df336941','14481388':'b2ffea8dce960b8aa1cbd86b88f8d160f2bb373562f4e9e2f6ecd13cca2d036e'}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main(parquet,image_csv,devdir,out):
 packet=guard();parquet=Path(parquet);image_csv=Path(image_csv);devdir=Path(devdir);out=Path(out)
 if sha(parquet)!=DATA_SHA or sha(image_csv)!=IMAGE_SHA:raise ValueError('source mismatch')
 if out.exists():raise FileExistsError(out)
 baseline={}
 for sid,digest in DEV_SHA.items():
  f=devdir/f'{sid}-development.json'
  if sha(f)!=digest:raise ValueError('baseline mismatch '+sid)
  baseline[sid]=json.loads(f.read_text(encoding='utf-8'))['result']
 control_rows=list(csv.DictReader((ROOT.parent/'vast-tde-expanded/control-metadata.csv').open(newline='',encoding='utf-8')))
 meta={}
 with image_csv.open(newline='',encoding='utf-8') as f:
  for x in csv.DictReader(f):meta[int(x['image_id'])]=(x['beam_bmaj'],x['beam_bmin'])
 out.mkdir(parents=True);manifest={'source_parquet_sha256':DATA_SHA,'image_csv_sha256':IMAGE_SHA,'development_sha256':DEV_SHA,'projections':{}}
 for sid,conf in packet['sources'].items():
  for phase in ('development','holdout'):
   groups=conf[phase];ids={i for group in groups for i in group}
   if phase=='development':controls={int(x['control_source']) for x in control_rows if x['target_source']==sid}
   else:controls={int(k) for k in baseline[sid]['control_baselines']}
   sources=controls|{int(sid)}
   tab=pq.read_table(parquet,columns=COLUMNS,filters=[('source','in',list(sources)),('image_id','in',list(ids))])
   if not set(tab.column('source').to_pylist())<=sources or not set(tab.column('image_id').to_pylist())<=ids:raise ValueError('unauthorized returned row')
   rows=tab.to_pylist()
   for r in rows:r['beam_bmaj'],r['beam_bmin']=meta[int(r['image_id'])]
   phase_rows(rows,ids,sources) # checks duplicate and membership; raw failed rows remain in projection
   name=f'{sid}-{phase}.jsonl.gz';path=out/name
   with path.open('xb') as f:
    with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:
     for r in rows:z.write((json.dumps(r,sort_keys=True,allow_nan=True)+'\n').encode())
   manifest['projections'][name]={'rows':len(rows),'sha256':sha(path),'source_ids':len(sources),'image_ids':sorted(ids)}
 (out/'projection-manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--parquet',required=True);ap.add_argument('--images',required=True);ap.add_argument('--development-dir',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();main(a.parquet,a.images,a.development_dir,a.output)
