"""Full synthetic geometry-to-selection injection; never accesses real FITS."""
import numpy as np
from astropy.io import fits
from scipy.ndimage import binary_dilation
from footprints import footprints
from development import assess
from prediction import construct
from validate_repeat import evaluate

h=fits.Header()
for k,v in {'DATE-BEG':'2026-08-13T16:25:09','DATE-AVG':'2026-08-13T16:27:22','DATE-END':'2026-08-13T16:29:37','TIMESYS':'UTC','BUNIT':'Gauss','DSUN_OBS':4.75e10,'DSUN_AU':.3175,'HGLN_OBS':-89.48,'HGLT_OBS':-15.75,'RSUN_REF':695700000.,'CTYPE1':'HPLN-TAN','CTYPE2':'HPLT-TAN','CUNIT1':'arcsec','CUNIT2':'arcsec','CRPIX1':128.,'CRPIX2':128.,'CRVAL1':0.,'CRVAL2':0.,'CDELT1':30.,'CDELT2':30.,'PXBEG1':65,'PXBEG2':65,'NAXIS1':256,'NAXIS2':256,'FILENAME':'synthetic'}.items():h[k]=v
patches=footprints(h,[h,h,h,h]);assert len(patches)==25 and all(x['eligible'] for x in patches)
first=np.zeros((256,256));second=first.copy()
blank=assess(first,second,patches);assert blank['status']=='no_candidate'
p=patches[12];m=p['maps'][1];mask=np.zeros_like(second,dtype=bool)
mask[np.floor(m['y']).astype(int),np.floor(m['x']).astype(int)]=True
second[binary_dilation(mask,iterations=4)]=20.
result=assess(first,second,patches);assert result['status']=='candidate_proposed' and result['candidate']==p['patch_id']
# Prediction consumes only a development report and never accesses FITS.
result.update({'holdout_pixels_read':False,'run_id':'toy','code_sha':'synthetic','packet_sha256':'synthetic','development_obs_ids':['early','late'],'holdout_obs_ids':['a','b'],'holdout_sha256':['a','b'],'tracked_footprint':[{'filename':n,**{k:np.asarray(v).tolist() for k,v in m.items()}} for n,m in zip(['early','late','a','b'],p['maps'])],'detector_frame_control':{'baseline_median':0.,'maps':[{'filename':n,'eligible':True,**{k:np.asarray(v).tolist() for k,v in m.items()}} for n,m in zip(['early','late','a','b'],p['maps'])]}})
pr=construct(result,'synthetic');assert pr['predicted_sign']==1 and pr['patch_id']==p['patch_id']
assert evaluate(pr,[second,second])['both_pass']
assert not evaluate(pr,[second,first])['both_pass']
print('synthetic 25-patch WCS mapping, null, extended injection, signed prediction and held-out falsification pass')
