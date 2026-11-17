"""Post-screen adapter: project authorized MMS1 records to a compact replay input."""
from pathlib import Path
import sys,json,hashlib
import cdflib,numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'frozen'))
from development import read_sc1,EXPECTED

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 if len(sys.argv)!=4:raise SystemExit('usage: make_projection.py MMS1_FGM_CDF MMS1_DIS_CDF NEW_OUTPUT_NPZ')
 f,d,out=map(Path,sys.argv[1:]); expected=json.loads((HERE/'evidence/independent.json').read_text())['source_hashes']
 for p in (f,d):
  if p.name not in expected or sha(p)!=expected[p.name]:raise ValueError('source hash mismatch '+p.name)
 if out.exists():raise FileExistsError(out)
 coverage={};origin,t,b,v,e,good=read_sc1(cdflib.CDF(str(f)),cdflib.CDF(str(d)),EXPECTED,coverage)
 with out.open('xb') as h:np.savez_compressed(h,origin_tt2000=origin,t=t,field=b,good=good)
 print(json.dumps({'sha256':sha(out),'source_hashes':expected,'coverage':coverage,'columns':['origin_tt2000','t_seconds','field_GSE_nT','aligned_quality_valid'], 'scope':'post-screen projection from authorized MMS1 only'},indent=2))
if __name__=='__main__':main()
