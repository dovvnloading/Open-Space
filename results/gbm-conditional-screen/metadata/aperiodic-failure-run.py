import numpy as np,json,sys,hashlib
from pathlib import Path
sys.path.insert(0,'D:/Open-Space/worktrees/gbm-conditional-qpo-20261117/results/gbm-conditional-screen')
from envelope import joint_bootstrap
rng=np.random.default_rng(909)
def burst(length):
 t=np.arange(length)/64;centers=[.35,.82,1.37,2.11,2.75,3.31,4.08,5.21,6.8,8.3,10.2,12.7,15.9,18.4];widths=[.045,.08,.03,.11,.05,.075,.08,.1,.04,.11,.075,.05,.09,.06];m=np.ones(length)*8
 for c,w in zip(centers,widths):m+=40*np.exp(-.5*((t-c)/w)**2)
 return np.stack([rng.poisson(m),rng.poisson(.7*m)])
x,y=burst(512),burst(1280);o={}
for width in (9,33,65):
 a=joint_bootstrap(np.random.default_rng(1000+width),[x,y],19,width);o[str(width)]={'p':a['p'],'observed':a['observed'],'n':19,'scope':a['scope']}
p=Path('D:/Open-Space/campaigns/gbm-envelope-synthetic-20261117/result.json');p.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2));print('sha',hashlib.sha256(p.read_bytes()).hexdigest())
