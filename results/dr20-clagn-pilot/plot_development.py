"""Development-only diagnostic plot for the nominal model outlier."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
src=Path(sys.argv[1]);rec=json.loads(Path(sys.argv[2]).read_text());sid=62239511;items=[r for r in rec if r['sdss_id']==sid and r['kind']!='hold'];z=next(r['z'] for r in items if r['kind']=='baseline');fig,ax=plt.subplots(2,1,figsize=(10,5),sharex=True)
for i,kind in enumerate(('baseline','dev')):
 r=next(r for r in items if r['kind']==kind);p=src/r['filename'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
 with fits.open(p) as h:
  c=h['COADD'].data;w=10**c['loglam']/(1+z);m=(w>=4750)&(w<=5100)&(c['ivar']>0)&(c['and_mask']==0);ax[i].plot(w[m],c['flux'][m],lw=.7);ax[i].set_ylabel('flux');ax[i].set_title(kind+' MJD '+str(r['mjd_list'][0]));ax[i].axvline(4862.68,color='gray',ls=':');ax[i].axvline(5008.239,color='gray',ls=':')
ax[-1].set_xlabel('Baseline-redshift rest wavelength (Å)');fig.tight_layout();fig.savefig(Path(__file__).with_name('outlier-development.png'),dpi=140)
