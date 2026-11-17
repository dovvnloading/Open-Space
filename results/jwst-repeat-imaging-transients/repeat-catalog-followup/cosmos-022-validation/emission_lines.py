"""Fit narrow optical lines in the pre-event MUSE nuclear aperture spectrum."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import least_squares
from scipy.special import ndtr

D=Path(__file__).resolve().parent
s=pd.read_csv(D/'muse_spectrum.csv.gz');wave=s.wavelength_A.values;flux=s.flux.values
errors=s.formal_error.values;dl=np.median(np.diff(wave));z=.5067
# Air wavelengths for MUSE AWAV.
lines={'NeV3426':[3425.88],'OII3727':[3726.03,3728.82],'NeIII3869':[3868.76],
       'Hdelta':[4101.74],'Hgamma':[4340.47],'HeII4686':[4685.68],
       'Hbeta':[4861.33],'OIII4959':[4958.91],'OIII5007':[5006.84]}
def fit(name,rest):
    centers=np.array(rest)*(1+z);mid=centers.mean();q=(abs(wave-mid)<35)&np.isfinite(flux)&np.isfinite(errors)&(errors>0)
    x=wave[q];y=flux[q];e=errors[q];side=abs(x-mid)>12
    def design(dx,sigma):
        cols=[np.ones(len(x)),x-mid]
        for c in centers:cols.append((ndtr((x+dl/2-c-dx)/sigma)-ndtr((x-dl/2-c-dx)/sigma))/dl)
        return np.column_stack(cols)
    def solve(dx,sigma):
        mat=design(dx,sigma);coef=np.linalg.lstsq(mat/e[:,None],y/e,rcond=None)[0]
        return mat,coef
    def residual(v):
        mat,co=solve(*v);return (mat@co-y)/e
    # Strong OII and OIII establish a common shift/width. The fainter lines use
    # a fixed template to avoid upward bias when testing their presence.
    if name in ('OII3727','OIII5007'):
        zfit=least_squares(residual,[0.,1.3],bounds=([-2.5,.8],[2.5,4])).x
    else:zfit=np.array([0.,1.5])
    m,c=solve(*zfit)
    delta=(y-m@c)/e
    scale=max(1.,1.4826*np.median(abs(delta[side]-np.median(delta[side]))))
    cov=np.linalg.inv((m/e[:,None]).T@(m/e[:,None]))*scale**2
    f=c[2:].sum();df=cov[2:,2:].sum()**.5
    return dict(line=name,flux_1e20=f,error_1e20=df,snr=f/df,delta_A=zfit[0],sigma_A=zfit[1],continuum_1e20_per_A=c[0],empirical_error_scale=scale)
if __name__=='__main__':
 rows=[fit(k,v) for k,v in lines.items()];out=pd.DataFrame(rows);out.to_csv(D/'muse_lines.csv',index=False)
 print(out[['line','flux_1e20','error_1e20','snr','delta_A','sigma_A']].round(2).to_string(index=False))
