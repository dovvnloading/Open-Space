"""Free the [O III] doublet ratio and fit the three-line common-redshift hypothesis to the grade-A spectrum."""
import json
import sys
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from infer import HERE

def run(cache):
    data=ET.parse(cache/'2655306065654633189.vot')
    v=np.array([[float(c.text or 'nan')for c in row]for row in data.findall('.//{*}TR')])
    w,f,mask,quality,var,ndith,err=v.T
    z0=2.540083647
    rest=np.array([4862.68,4960.295,5008.24])
    centres=rest*(1+z0)
    keep=(mask.astype(int)%2==0)&(mask<64)&(ndith>=2)&(err>0)&np.isfinite(f+err)
    keep&=((abs(w-centres[0])<165)|((w>centres[1]-115)&(w<centres[2]+180)))
    w,f,e=w[keep],f[keep],err[keep]
    blue=w<(centres[0]+centres[1])/2
    def predicted(p):
        z,sigma=p[:2]
        fluxes=p[2:5]
        cont=np.where(blue,p[5]+p[6]*(w-centres[0])/160,p[7]+p[8]*(w-centres[2])/260)
        templates=np.exp(-.5*((w[:,None]-rest*(1+z))/sigma)**2)/(sigma*np.sqrt(2*np.pi))
        return cont+templates@fluxes
    fit=least_squares(lambda p:(predicted(p)-f)/e,[z0,28,2,3,9,.05,0,.05,0],
                      bounds=([z0-.015,8,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf],
                              [z0+.015,90,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))
    p=fit.x
    chi2=float(fit.fun@fit.fun)
    cov=np.linalg.inv(fit.jac.T@fit.jac)*max(1,chi2/(len(w)-len(p)))
    sigma=np.sqrt(np.diag(cov))
    ratio=p[4]/p[3]
    j=np.zeros_like(p);j[4]=1/p[3];j[3]=-ratio/p[3]
    result=dict(object_id='2655306065654633189',n_bins=len(w),n_fit=len(p),chi2=chi2,z_fit=p[0],z_formal_error=sigma[0],
                FWHM_A=2.35482*p[1],flux_unit='1e-16 erg s^-1 cm^-2',
                line_names=['Hbeta','[OIII]4960','[OIII]5008'],line_fluxes=p[2:5].tolist(),line_errors=sigma[2:5].tolist(),
                SNR=(p[2:5]/sigma[2:5]).tolist(),OIII_5008_to_4960=ratio,ratio_error=float(np.sqrt(j@cov@j)))
    (HERE/'sir_gradeA_features.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    run(Path(sys.argv[1]))
