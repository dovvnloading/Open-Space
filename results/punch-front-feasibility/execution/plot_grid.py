import json,pathlib,hashlib,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=pathlib.Path('D:/Open-Space/campaigns/punch-development-16f73d35-20261117/development.json'); x=json.loads(p.read_text())['result']; g=np.array([[np.nan if v is None else v for v in row] for row in x['difference_grid_64']]); q=np.nanpercentile(abs(g),98)
fig,ax=plt.subplots(figsize=(7,6)); z=ax.imshow(g,origin='lower',cmap='RdBu_r',vmin=-q,vmax=q,interpolation='nearest'); fig.colorbar(z,ax=ax,label='Block-mean difference (header brightness units)'); ax.set(title='PUNCH 00:08 − 00:00 development difference\n64×64 descriptive grid; uncalibrated differences',xlabel='Frame-0 block column',ylabel='Frame-0 block row'); fig.tight_layout(); fig.savefig('D:/Open-Space/campaigns/punch-grid-review-16f73d35-20261117/grid.png',dpi=160); print('input',hashlib.sha256(p.read_bytes()).hexdigest(),'finite_blocks',int(np.isfinite(g).sum()),'scale98',q)
