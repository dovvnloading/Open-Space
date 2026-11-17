"""Combine two frozen replay stages into an exclusive external directory."""
import argparse,shutil,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from replay_guard import absolute_external,fresh_dirs,verify_seal,seal

def main():
 p=argparse.ArgumentParser();p.add_argument('--screen',required=True,type=Path);p.add_argument('--calibration',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
 inputs=[absolute_external(a.screen),absolute_external(a.calibration)];dest=absolute_external(a.out)
 if any(x==dest or x in dest.parents or dest in x.parents for x in inputs):raise ValueError('input/output overlap')
 names=[verify_seal(inputs[0],'screen'),verify_seal(inputs[1],'calibration')]
 if set(names[0])&set(names[1]):raise ValueError('stage filenames collide')
 (dest,)=fresh_dirs(dest)
 for src,files in zip(inputs,names):
  for name in files:shutil.copyfile(src/name,dest/name)
 seal(dest,'combined')
if __name__=='__main__':main()
