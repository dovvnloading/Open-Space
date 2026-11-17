import numpy as np,json,sys,hashlib
from pathlib import Path
sys.path.insert(0,'D:/Open-Space/worktrees/gbm-conditional-qpo-20261117/results/gbm-conditional-screen');from envelope import joint_bootstrap
rng=np.random.default_rng(5601);out=[]
for case in range(20):
 def burst(length):
  t=np.arange(length)/64;mean=np.ones(length)*8
  for _ in range(12):
   center=rng.uniform(0,t[-1]);width=rng.uniform(.025,.16);amp=rng.uniform(12,65);mean+=amp*np.exp(-.5*((t-center)/width)**2)
  return np.stack([rng.poisson(mean),rng.poisson(mean*.7)])
 x,y=burst(512),burst(1280);r=joint_bootstrap(np.random.default_rng(7000+case),[x,y],19,33);out.append({'case':case,'p':r['p'],'exceedances':int(round(r['p']*20-1))})
p=Path('D:/Open-Space/campaigns/gbm-pulse-ensemble-20261117/result.json');p.write_text(json.dumps({'n_cases':20,'simulations_per_case':19,'width':33,'cases':out,'scope':'synthetic irregular pulse population; Monte Carlo floor .05, not real data calibration'},indent=2)+'\n');print('min_p_cases',sum(z['p']==.05 for z in out),'of',len(out),'sha',hashlib.sha256(p.read_bytes()).hexdigest())
