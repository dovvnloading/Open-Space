"""Synthetic-only selection, registration and missing-data controls."""
import numpy as np
from development import assess,sample,read_dev

rng=np.random.default_rng(2901);base=rng.normal(0,.2,(128,128));repeat=base.copy();patches=[]
for i in range(5):
    for j in range(5):
        x=10+i*22;y=10+j*22
        xx,yy=np.meshgrid(np.linspace(x-2,x+2,5),np.linspace(y-2,y+2,5))
        m={'x':xx.ravel(),'y':yy.ravel(),'mu':np.ones(25)}
        patches.append({'patch_id':f'g{i}{j}','eligible':True,'reasons':[],'maps':[m,m]})
zero=assess(base,repeat,patches);assert zero['status']=='no_candidate'
# A repeatable extended enhancement survives all +/-1-pixel controls.
repeat[72:81,50:59]+=20.
injected=assess(base,repeat,patches);assert injected['status']=='candidate_proposed' and injected['candidate']=='g23'
# A one-pixel defect fails the registration robustness gate.
repeat=base.copy();repeat[76,54]+=500.
spike=assess(base,repeat,patches);assert spike['status']=='no_candidate'
# A NaN in the bilinear footprint must exclude its patch, not be interpreted as zero.
repeat=base.copy();repeat[10,10]=np.nan
bad=assess(base,repeat,patches);assert bad['invalid_stencils']==1 and bad['rows'][0]['status']=='invalid_stencil'
print('synthetic flat, extended injection, registration spike and invalid-stencil controls pass')

# Simulated independent Gaussian noise, not a calibration of real solar data.
# Predeclared trial count; report measured outcomes without significance claim.
null_hits=0;injection_hits=0
for _ in range(40):
    first=rng.normal(0,.5,(128,128));second=rng.normal(0,.5,(128,128))
    null_hits+=assess(first,second,patches)['status']=='candidate_proposed'
    second[72:81,50:59]+=20.
    injected=assess(first,second,patches)
    injection_hits+=(injected['candidate']=='g23')
print('Gaussian toy null hits',null_hits,'/40; extended injection recovery',injection_hits,'/40')
assert null_hits==0 and injection_hits>=38

try:read_dev(('solo_L2_phi-fdt-blos_20260815T162509_V03_0648150506.fits','unused'))
except ValueError as exc:assert 'development filenames' in str(exc)
else:raise AssertionError('holdout filename accepted by development reader')
print('holdout filename rejected before file access')
