"""Read-only preflight for a future retrospective pixel runner; no image decoding."""
import argparse, hashlib, json
from pathlib import Path

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 return h.hexdigest()

def verify_mapping(path, expected):
 if digest(path)!=expected: raise ValueError('mapping digest mismatch')
 data=json.loads(path.read_text(encoding='utf-8'))
 rows=data['rows'];ids={(r['experiment'],r['draw']) for r in rows}
 if len(rows)!=800 or len(ids)!=800:raise ValueError('invalid draw membership')
 if any(r['status']!='unique' or len(r['matches'])!=1 for r in rows):raise ValueError('unmapped or ambiguous draw; manual disposition required')
 if {e:sum(r['experiment']==e for r in rows) for e in ('published_centroids','centroid_stress')}!={'published_centroids':400,'centroid_stress':400}:raise ValueError('wrong experiment counts')
 return data

def preflight(scenario,mapping,mapping_sha,inputs,output):
 if scenario not in ('historical','dec1','dec2'):raise ValueError('invalid scenario')
 if output.exists():raise FileExistsError(output)
 data=verify_mapping(mapping,mapping_sha)
 manifest=json.loads(inputs.read_text(encoding='utf-8'))
 if not manifest.get('inputs'):raise ValueError('empty input manifest')
 for entry in manifest['inputs']:
  path=Path(entry['path'])
  if not path.is_absolute() or not path.is_file() or path.stat().st_size!=entry['size'] or digest(path)!=entry['sha256']:raise ValueError('input verification failed: '+str(path))
 # Actual runners are intentionally not connected until reviewed.
 return {'scenario':scenario,'verified_rows':len(data['rows']),'verified_inputs':len(manifest['inputs']),'status':'preflight only; pixel runner not implemented'}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--scenario',required=True,choices=['historical','dec1','dec2']);p.add_argument('--mapping',required=True,type=Path);p.add_argument('--mapping-sha256',required=True);p.add_argument('--inputs',required=True,type=Path);p.add_argument('--cache',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
 if a.cache.exists() and not a.cache.is_dir():raise ValueError('cache must be a directory')
 print(json.dumps(preflight(a.scenario,a.mapping,a.mapping_sha256,a.inputs,a.out),indent=2))
if __name__=='__main__':main()
