"""Freeze/verify the DES timing-case executable packet before image decoding."""
import argparse, hashlib, json, platform, subprocess, sys
from pathlib import Path
import astropy, erfa, numpy, scipy
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
FILES=['archives/decam_difference.py','archives/planet_injections.py','archives/decam_inputs.json','orbits/orbit_input.json','orbits/propagate.py','orbits/timing_scenarios.py','timing_replay/des_runner.py','timing_replay/preflight.py','timing_replay/scenario_states.py','timing_replay/controls.py','timing_replay/des_packet.py']
def digest(b):return hashlib.sha256(b).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=REPO)
def runtime():return {'python':sys.version,'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'astropy':astropy.__version__,'erfa':erfa.__version__}
def snapshot(phase='historical'):
 if phase not in ('historical','dec1','dec2'):raise ValueError('invalid phase')
 head=git('rev-parse','HEAD').decode().strip();files={}
 for rel in FILES:
  path=ROOT/rel;blob=git('show',head+':'+path.relative_to(REPO).as_posix());disk=path.read_bytes()
  if blob!=disk:raise ValueError('worktree differs from Git blob: '+rel)
  files[rel]=digest(blob)
 return {'schema':1,'phase':phase,'git_head':head,'files':files,'runtime':runtime()}
def verify(manifest_path,manifest_sha,auth_path,auth_sha,scenario,mapping,mapping_sha,cache,out):
 if digest(manifest_path.read_bytes())!=manifest_sha or digest(auth_path.read_bytes())!=auth_sha:raise ValueError('packet or authorization hash mismatch')
 manifest=json.loads(manifest_path.read_bytes());auth=json.loads(auth_path.read_bytes())
 if manifest!=snapshot(scenario):raise ValueError('code, runtime or Git head differs from frozen packet')
 expected={'phase':scenario,'manifest_sha256':manifest_sha,'mapping_sha256':mapping_sha,'mapping':str(mapping.resolve()),'cache':str(cache.resolve()),'out':str(out.resolve())}
 if scenario not in ('historical','dec1','dec2') or auth!=expected:raise ValueError('unauthorized phase or configuration')
 return manifest

def main():
 a=argparse.ArgumentParser();a.add_argument('--freeze',required=True,type=Path);a.add_argument('--scenario',required=True,choices=['historical','dec1','dec2']);q=a.parse_args()
 if q.freeze.exists():raise FileExistsError(q.freeze)
 q.freeze.parent.mkdir(parents=True,exist_ok=True)
 with q.freeze.open('x',encoding='utf-8',newline='\n') as f:json.dump(snapshot(q.scenario),f,indent=2);f.write('\n')
 print(digest(q.freeze.read_bytes()))
if __name__=='__main__':main()
