"""Bind tracked packet, four exact downloaded TTE bytes and runtime; no EVENTS decoding."""
from pathlib import Path
import subprocess,hashlib,json,platform,sys
import numpy,astropy
P=Path(__file__).resolve().parent;ROOT=P.parents[1];RAW=Path('D:/Open-Space/campaigns/gbm-tte-metadata-20261117/raw');RECEIPTS=RAW.parent/'receipts.json'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def snapshot():
 head=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip();names=subprocess.check_output(['git','-C',str(ROOT),'ls-files','results/gbm-conditional-screen'],text=True).splitlines();items={}
 for name in names:
  p=ROOT/name;raw=p.read_bytes()
  if raw!=subprocess.check_output(['git','-C',str(ROOT),'show',head+':'+name]):raise ValueError('packet differs from git '+name)
  items[name]=hashlib.sha256(raw).hexdigest()
 receipts=json.loads(RECEIPTS.read_text())['files'];files={}
 for x in receipts:
  name=x['url'].split('/')[-1];p=RAW/name
  if sha(p)!=x['sha256']:raise ValueError('input differs '+name)
  files[name]=x['sha256']
 if len(files)!=4:raise ValueError('expected 4 TTEs')
 return {'git_head':head,'packet_sha256':items,'tte_sha256':files,'receipts_sha256':sha(RECEIPTS),'python':platform.python_version(),'numpy':numpy.__version__,'astropy':astropy.__version__,'statistic':'joint maximum; no separate independent confirmation'}
def verify(path,expected):
 if sha(path)!=expected or json.loads(Path(path).read_text())!=snapshot():raise ValueError('manifest mismatch')
 return json.loads(Path(path).read_text())
if __name__=='__main__':
 p=Path(sys.argv[1]);raw=(json.dumps(snapshot(),indent=2,sort_keys=True)+'\n').encode()
 with p.open('xb') as h:h.write(raw)
 print(hashlib.sha256(raw).hexdigest())
