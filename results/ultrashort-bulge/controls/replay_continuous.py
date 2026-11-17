"""Replay only the known control's continuous fits, without native grid execution."""
import argparse,gzip,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import earth_control
from replay_guard import external,fresh_for_inputs

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 source=external(a.input);raw=source.read_bytes()
 if raw[:2]==b'\x1f\x8b':raw=gzip.decompress(raw)
 if hashlib.sha256(raw).hexdigest()!='819169bc4c6ab82941c2c01b48c8ba9c331bd55735348d30f2d5715bc2cb8c0c':raise ValueError('control input mismatch')
 out=fresh_for_inputs(a.out,source);setup,screen=earth_control.fit.setup,earth_control.fit.screen
 try:
  earth_control.fit.setup=lambda *args:None
  earth_control.fit.screen=lambda *args,**kw:{'status':'not replayed; C++ grid bypassed'}
  result=earth_control.run(raw)
 finally:
  earth_control.fit.setup=setup;earth_control.fit.screen=screen
 (out/'continuous.json').write_text(json.dumps(result,indent=2,default=lambda x:x.item())+'\n')
if __name__=='__main__':main()
