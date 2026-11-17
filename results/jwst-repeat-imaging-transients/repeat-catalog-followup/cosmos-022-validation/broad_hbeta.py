"""Test historical broad H-beta at two fixed widths, with local linear continuum."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import ndtr
D=Path(__file__).resolve().parent
t=pd.read_csv(D/'muse_spectrum.csv.gz')
c=4861.33*1.5067
sel=abs(t.wavelength_A-c)<120
x,y,e=(t.loc[sel,k].values for k in ['wavelength_A','flux','formal_error'])
g=np.isfinite(y)&np.isfinite(e)&(e>0)
x,y,e=x[g],y[g],e[g];step=np.median(np.diff(x))
def line(sig):
 return (ndtr((x+step/2-c)/sig)-ndtr((x-step/2-c)/sig))/step
rows=[]
for fwhm in (2000,4000):
 sig=c*(fwhm/299792.458)/2.35482
 a=np.column_stack([np.ones_like(x),x-c,line(1.50),line(sig)])
 fit=np.linalg.lstsq(a/e[:,None],y/e,rcond=None)[0]
 resid=(y-a@fit)/e
 inflate=max(1.,1.4826*np.median(abs(resid-np.median(resid))))
 cov=np.linalg.inv((a/e[:,None]).T@(a/e[:,None]))*inflate**2
 err=cov[-1,-1]**.5;f=fit[-1]
 rows.append(dict(fwhm_kms=fwhm,narrow_flux_1e20=fit[-2],broad_flux_1e20=f,broad_error_1e20=err,snr=f/err,
                 conservative_3sigma_upper_1e20=max(0,f)+3*err))
r=pd.DataFrame(rows);r.to_csv(D/'muse_broad_hbeta.csv',index=False);print(r.round(1).to_string(index=False))
