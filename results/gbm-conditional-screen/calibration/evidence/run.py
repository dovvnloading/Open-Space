import sys,json,hashlib
from pathlib import Path
sys.path.insert(0,'D:/Open-Space/worktrees/gbm-synthetic-null-20261117/results/gbm-retrospective-null')
from calibration import summarize
out=summarize(100);p=Path('D:/Open-Space/campaigns/gbm-retrospective-null-20261117/result.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
