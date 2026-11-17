"""Render only previously authorized diagnostic stamps and empirical PSFs."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
with np.load(HERE/'evidence/diagnostic/run/native_stamps.npz') as z:
 fig,ax=plt.subplots(2,5,figsize=(12,5),constrained_layout=True)
 for row,case in enumerate(('dec1','dec2')):
  for col,(name,label) in enumerate((('shifted_epoch1','First exposure'),('shifted_epoch2','Shifted second'),('detector_fixed_epoch2','Detector-fixed second'),('psf1','First PSF'),('psf2','Second PSF'))):
   a=z[f'{case}_{name}'];v=np.nanpercentile(a,[5,95]);ax[row,col].imshow(a,origin='lower',cmap='gray',vmin=v[0],vmax=v[1]);ax[row,col].set_title(f'{case}: {label}',fontsize=9);ax[row,col].set_xticks([]);ax[row,col].set_yticks([])
 fig.suptitle('Fixed DES positions and empirical PSFs; same exposures, per-panel contrast scaling',fontsize=10)
 fig.savefig(HERE/'fixed-feature.png',dpi=150);plt.close(fig)
