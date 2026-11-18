"""Run the deposited search in a fresh D: campaign with verified input bytes."""
import argparse,hashlib,importlib.metadata,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def verify(root,manifest):
 for name,expected in manifest['files'].items():
  p=(root/name).resolve()
  if not p.is_relative_to(root.resolve()) or not p.is_file() or sha(p)!=expected: raise ValueError('input/code mismatch: '+name)
 actual={'python':sys.version.split()[0],**{x:importlib.metadata.version(x) for x in ('numpy','pandas','scipy')}}
 if actual!=manifest['runtime']: raise ValueError('runtime mismatch: '+str(actual))
 return actual

def create(root,output,manifest,stage):
 output=Path(output)
 if not output.is_absolute(): raise ValueError('absolute output required')
 resolved=output.resolve()
 if not (str(resolved).lower().startswith('d:\\') or str(resolved).startswith('/mnt/d/')): raise ValueError('output must resolve under D:')
 if resolved.is_relative_to(root.resolve()): raise ValueError('source tree output forbidden')
 verify(root,manifest)
 resolved.mkdir(parents=False,exist_ok=False)
 for name in manifest['files']:
  src=root/name; dest=resolved/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
 (resolved/'measurement').mkdir()
 env=os.environ.copy();env['HALO_REPLAY_OUTPUT']=str(resolved/'measurement');env['HALO_REPLAY_ISOLATED']='1'
 stages=[['search.py']] if stage=='search' else [['search.py'],['off_control.py'],['inject.py'],['inject.py','--off-control']]
 results=[]
 for args in stages:
  r=subprocess.run([sys.executable,str(resolved/'scripts'/args[0]),*args[1:]],cwd=resolved,env=env)
  results.append({'args':args,'exit_code':r.returncode})
  if r.returncode: break
 (resolved/'replay-receipt.json').write_text(json.dumps({'stage':stage,'manifest_sha256':sha(root/'replay_manifest.json'),'steps':results},indent=2)+'\n')
 if any(x['exit_code'] for x in results): raise RuntimeError('replay failed; inspect receipt')

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--stage',choices=('search','full'),required=True);a=p.parse_args()
 if sha(ROOT/'replay_manifest.json')!=a.manifest_sha256: raise ValueError('manifest digest mismatch')
 m=json.loads((ROOT/'replay_manifest.json').read_text());create(ROOT,a.output,m,a.stage)
if __name__=='__main__':main()
