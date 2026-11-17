"""Post-screen adapter; read only the frozen MAG 00–12 UTC interval."""
import sys,json,hashlib
from pathlib import Path
import cdflib,numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'frozen'))
from development import mag_arrays

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 if len(sys.argv)!=3:raise SystemExit('usage: make_projection.py MAG_CDF NEW_OUTPUT.npz')
 p,out=map(Path,sys.argv[1:]);expected=json.loads((HERE/'evidence/independent.json').read_text())['source_sha256']
 if sha(p)!=expected or out.exists():raise ValueError('source mismatch or output exists')
 t,b,good,coverage=mag_arrays(cdflib.CDF(p))
 with out.open('xb') as f:np.savez_compressed(f,t=t,b=b,good=good)
 print(json.dumps({'source_sha256':expected,'sha256':sha(out),'coverage':coverage,'scope':'post-screen MAG first twelve UTC hours only'},indent=2))
if __name__=='__main__':main()
