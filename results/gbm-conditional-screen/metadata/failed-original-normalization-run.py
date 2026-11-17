import numpy as np,json,hashlib,sys
from pathlib import Path
sys.path.insert(0,'D:/Open-Space/worktrees/gbm-conditional-qpo-20261117/results/gbm-conditional-screen')
from simulator import monte_carlo
n=100;train=monte_carlo(417,n);test=monte_carlo(519,n);alt=monte_carlo(617,n,kind='pulse_train');injected=monte_carlo(719,n,injection=.15);cut=float(np.quantile(train,.95,method='higher'))
o={'n_each':n,'threshold_training_95pct':cut,'independent_null_crossings':int((test>=cut).sum()),'pulse_train_crossings':int((alt>=cut).sum()),'injected_crossings':int((injected>=cut).sum()),'scope':'synthetic toy calibration only; empirical rate uncertainty and model misspecification remain; not a real-burst tail probability'};p=Path('D:/Open-Space/campaigns/gbm-synthetic-20261117/synthetic-validation.json');p.write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
