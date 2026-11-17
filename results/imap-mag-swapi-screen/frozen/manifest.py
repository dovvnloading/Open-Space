"""Review packet metadata and SHA guards; no physical data values."""
import sys,platform,json,hashlib,subprocess
from pathlib import Path
import cdflib,numpy
P=Path(__file__).resolve().parent;ROOT=P.parents[1];RAW=Path('D:/Open-Space/campaigns/heliophysics-schema-20261117/raw')
FILES=['imap_mag_l2_norm-gse_20260731_v002.0003.cdf','imap_swapi_l3a_proton-sw_20260731_v002.0006.cdf']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def snapshot():
 head=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip();packet={}
 for p in sorted(P.iterdir()):
  if p.is_file():
   rel=p.relative_to(ROOT).as_posix();assert p.read_bytes()==subprocess.check_output(['git','-C',str(ROOT),'show',head+':'+rel]),rel;packet[p.name]=sha(p)
 return {'git_head':head,'packet':packet,'input_hashes':{n:sha(RAW/n) for n in FILES},'python':platform.python_version(),'numpy':numpy.__version__,'cdflib':cdflib.__version__,'phase':'MAG 00-12 only; SWAPI sealed'}
def verify(p,digest):
 if sha(p)!=digest or json.loads(Path(p).read_text())!=snapshot():raise ValueError('manifest mismatch')
if __name__=='__main__':
 p=Path(sys.argv[1]);raw=(json.dumps(snapshot(),indent=2,sort_keys=True)+'\n').encode()
 with p.open('xb') as f:f.write(raw)
 print(hashlib.sha256(raw).hexdigest())
