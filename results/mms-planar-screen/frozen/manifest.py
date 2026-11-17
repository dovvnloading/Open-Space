"""Bind exact packet bytes, inputs, metadata and interpreter before development."""
import hashlib,json,platform,subprocess,sys
from pathlib import Path
import cdflib,numpy
P=Path(__file__).resolve().parent;ROOT=P.parents[1];RAW=Path('D:/Open-Space/campaigns/heliophysics-schema-20261117/raw');META=Path('D:/Open-Space/campaigns/mms-planar-metadata-20261117/metadata.json')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def snapshot():
 head=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip();files=sorted(p for p in P.iterdir() if p.is_file())
 d={}
 for p in files:
  rel=p.relative_to(ROOT).as_posix();blob=subprocess.check_output(['git','-C',str(ROOT),'show',head+':'+rel]);
  if blob!=p.read_bytes():raise ValueError('packet differs from Git blob: '+rel)
  d[p.name]=sha(p)
 records=json.loads(META.read_text())['records'];inputs={}
 for name,r in records.items():
  p=RAW/name
  if sha(p)!=r['sha256']:raise ValueError('input digest: '+name)
  inputs[name]=r['sha256']
 if len(inputs)!=8:raise ValueError('eight inputs required')
 return {'git_head':head,'packet':d,'inputs':inputs,'metadata_sha256':sha(META),'python':platform.python_version(),'numpy':numpy.__version__,'cdflib':cdflib.__version__,'source':'MMS1','heldout':['MMS2','MMS3','MMS4']}
def verify(path,expected):
 p=Path(path)
 if sha(p)!=expected or json.loads(p.read_text())!=snapshot():raise ValueError('manifest mismatch')
 return json.loads(p.read_text())
if __name__=='__main__':
 p=Path(sys.argv[1]);raw=(json.dumps(snapshot(),sort_keys=True,indent=2)+'\n').encode();
 with p.open('xb') as f:f.write(raw)
 print(hashlib.sha256(raw).hexdigest())
