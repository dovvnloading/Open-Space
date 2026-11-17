"""Replay deposited synthetic ensemble from fixed seeds; no observations accessed."""
import hashlib,json
from pathlib import Path
from calibration import summarize
P=Path(__file__).resolve().parent;root=P.parent
expected={'executed/calibration.py':'83b06b793a4d728c32496dfb6fe5533b29f542911de9ed7cffbc23bdb1c01449','../simulator.py':'04a864ee2555c50af7a2c2bca11f8294811675ab6654952476bb8ffa69d8e747','../envelope.py':'ec20828da8f803e6a7c5d29a01da0edfa1e0f6e495d0afc4a22c4cb05def0d69'}
for name,digest in expected.items():
 if hashlib.sha256((P/name).read_bytes()).hexdigest()!=digest:raise ValueError('code digest mismatch: '+name)
a=summarize(100);b=json.loads((P/'evidence/result.json').read_text())
if a!=b:raise ValueError('synthetic result mismatch')
print(json.dumps(a,indent=2))
