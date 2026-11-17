"""Freeze and authorize the fixed DES candidate diagnostic before FITS decode."""
import argparse,hashlib,json,platform,subprocess,sys
from pathlib import Path
import astropy,erfa,numpy,scipy
from feature_diagnostic import ROOT,LABELS,CASES,sha
REPO=ROOT.parents[1]
FILES=['archives/decam_difference.py','archives/decam_inputs.json','orbits/propagate.py','orbits/orbit_input.json','orbits/timing_scenarios.py','timing_replay/feature_diagnostic.py','timing_replay/feature_packet.py']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def snapshot(cache,prior,out):
 for p in (cache,prior,out):
  if not p.is_absolute() or p.resolve().drive.upper()!='D:':raise ValueError('D absolute path required')
 if out.exists():raise FileExistsError(out)
 head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip();code={}
 for rel in FILES:
  p=ROOT/rel;blob=subprocess.check_output(['git','show',head+':'+p.relative_to(REPO).as_posix()],cwd=REPO)
  if p.read_bytes()!=blob:raise ValueError('worktree/blob mismatch '+rel)
  code[rel]=hashlib.sha256(blob).hexdigest()
 receipts=json.loads((ROOT/'archives/decam_inputs.json').read_bytes());inputs={}
 for label in LABELS:
  matches=[r for r in receipts if r['label']==label]
  if len(matches)!=1:raise ValueError('receipt identity')
  r=matches[0];p=cache/(label+'.response')
  if not p.is_file() or p.stat().st_size!=r['size'] or sha(p)!=r['sha256']:raise ValueError('input mismatch '+label)
  inputs[label]={'sha256':r['sha256'],'bytes':r['size']}
 priors={}
 for case,c in CASES.items():
  d=prior/c['prior']/'run';priors[case]={x:digest(d/x) for x in ('difference_result.json','dense_peaks.json','difference_tracks.csv.gz','run_manifest.json','run_status.json')}
 return {'phase':'fixed_peak_diagnostic','git_head':head,'code_sha256':code,'inputs':inputs,'prior_sha256':priors,'cache':str(cache.resolve()),'prior_root':str(prior.resolve()),'out':str(out.resolve()),'runtime':{'python':sys.version,'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'astropy':astropy.__version__,'erfa':erfa.__version__}}
def verify(manifest,msha,auth,asha,cache,prior,out):
 if digest(manifest)!=msha or digest(auth)!=asha:raise ValueError('authorization or manifest digest')
 frozen=json.loads(manifest.read_bytes());approved=json.loads(auth.read_bytes())
 if frozen!=snapshot(cache,prior,out):raise ValueError('frozen packet mismatch')
 if approved!={'phase':'fixed_peak_diagnostic','manifest_sha256':msha,'out':str(out.resolve())}:raise ValueError('authorization mismatch')
def main():
 p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--prior-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--freeze',type=Path,required=True);a=p.parse_args()
 if a.freeze.exists():raise FileExistsError(a.freeze)
 record=snapshot(a.cache,a.prior_root,a.out);a.freeze.parent.mkdir(parents=True,exist_ok=True)
 with a.freeze.open('x',encoding='utf-8',newline='\n') as f:json.dump(record,f,indent=2);f.write('\n')
 print(digest(a.freeze))
if __name__=='__main__':main()
