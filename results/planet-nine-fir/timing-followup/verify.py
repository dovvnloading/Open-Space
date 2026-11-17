"""Offline integrity and bounded-count checks for the retrospective DES record."""
import collections,csv,gzip,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def data(path):return gzip.decompress(path.read_bytes()) if path.name.endswith('.json.gz') else path.read_bytes()
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 for line in (HERE/'SHA256SUMS').read_text().splitlines():
  digest,name=line.split('  ',1);assert sha((HERE/name).read_bytes())==digest,name
 mapping=json.loads(data(HERE/'evidence/mapping/mapping.json.gz'));assert len(mapping['rows'])==800
 assert collections.Counter(r['status'] for r in mapping['rows'])=={'unique':800}
 assert len({(r['experiment'],r['draw']) for r in mapping['rows']})==800
 for case in ('historical','dec1','dec2'):
  base=HERE/'evidence'/case;receipt=json.loads((base/'execution-receipt.json').read_bytes());outputs=receipt.get('outputs',receipt.get('run_output_sha256'))
  for name,digest in outputs.items():
   name=name.replace('\\','/').removeprefix('run/');p=base/'run'/name
   if not p.exists():p=p.with_name(p.name+'.gz')
   assert p.is_file(),str(p)
   assert sha(data(p))==digest,(case,name)
  status=json.loads((base/'run/run_status.json').read_bytes());assert status['status']=='complete'
  outcomes=json.loads(data(base/'run/draw_outcomes.json.gz'));assert len(outcomes)==800
  assert all(len(r['models'])==21 for r in outcomes)
  groups=collections.Counter((r['experiment'],r['state_status'],r['image_status']) for r in outcomes)
  expected=(389,11,379,21) if case=='historical' else (388,12,385,15)
  assert tuple(groups[(e,'bound',s)] for e,s in [('published_centroids','eligible'),('published_centroids','masked'),('centroid_stress','eligible'),('centroid_stress','masked')])==expected
  model=json.loads((base/'run/model_recovery.json').read_bytes()); dark={(r['experiment'],r['n_recovered']) for r in model if r['mass_Mearth']==5 and r['radius_Rearth']==2.71 and r['geometric_albedo_r']==0.1}
  assert dark==({'published_centroids':379,'centroid_stress':372}.items() if case=='historical' else {'published_centroids':375 if case=='dec1' else 374,'centroid_stress':374}.items())
 for case,code in [('historical','581c9c4e'),('dec1','f6b97d8d'),('dec2','f6b97d8d'),('diagnostic','ca3d18bf')]:
  packet=json.loads((HERE/'evidence'/f'{case}-packet'/'manifest.json').read_bytes())
  for name,digest in packet.get('files',packet.get('code_sha256',{})).items():
   assert sha((HERE/'executed-code'/code/Path(name.replace('\\','/'))).read_bytes())==digest,(case,name)
 diag=HERE/'evidence/diagnostic';receipt=json.loads((diag/'execution-receipt.json').read_bytes())
 for name,digest in receipt['output_sha256'].items():assert sha((diag/'run'/name).read_bytes())==digest,name
 d=json.loads((HERE/'evidence/diagnostic/run/diagnostic.json').read_bytes())
 for phase in ('dec1','dec2'):assert d[phase]['reconstruction']['status']=='passed'
 print('Frozen bytes, 800 identities, run outcomes and diagnostic reconstruction verified.')
if __name__=='__main__':main()
