"""Fit the original XSPEC line-emitting MEKAL table, with each elemental abundance at one solar."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
from fit_pha import read, sigma, profile_stat

MEKAL_URL='https://heasarc.gsfc.nasa.gov/FTP/software/xspec/spectral/modelData/mekal.mod'

def run(path,table,fixed_nh=None):
    e,m,n,b,alpha,_,_,rlo,rhi=read(path,.2)
    with fits.open(table) as f:
        s=f['SPECTRA'].data
        spec=np.array(s['INTPSPEC'],dtype=float)
        for k in s.columns.names:
            if k.startswith('ADDSP'):spec+=s[k]
        t=np.array(s['PARAMVAL'],dtype=float)
        if f['PARAMETERS'].data['METHOD'][0]!=0:raise ValueError('Expected linear kT interpolation')
        grid=np.r_[f['ENERGIES'].data['ENERG_LO'],f['ENERGIES'].data['ENERG_HI'][-1]].astype(float)
    # OGIP spectra are integrated in each native model bin. Redistribute
    # their photon *counts* to true RMF energy bins by cumulative sums.
    def regrid(source):
        cdf=np.r_[0.,np.cumsum(source)]
        return (np.interp(rhi,grid,cdf)-np.interp(rlo,grid,cdf))/(rhi-rlo)
    photons=np.maximum(np.array([regrid(row) for row in spec]),0.)
    absorb=sigma(e)*1e21
    def evaluate(x):
        total,nh,temp=np.exp(x[0]),x[1],np.exp(x[2])
        i=np.clip(np.searchsorted(t,temp)-1,0,len(t)-2)
        w=(temp-t[i])/(t[i+1]-t[i])
        unit=m@(np.exp(-nh*absorb)*(photons[i]+w*(photons[i+1]-photons[i])))
        return profile_stat(n,b,alpha,total*unit/max(unit.sum(),1e-200))
    limits=[(-8,16),(0,10) if fixed_nh is None else (fixed_nh,fixed_nh),(np.log(.02),np.log(5))]
    starts=np.geomspace(.025,3.5,9)
    fits_=[minimize(evaluate,[np.log(max(1,n.sum()-alpha*b.sum())),nh,np.log(kt)],bounds=limits,method='L-BFGS-B') for kt in starts for nh in [0.,1.,3.]]
    best=min(fits_,key=lambda x:x.fun)
    return {'file':path.name,'table':MEKAL_URL,'table_sha256':hashlib.sha256(table.read_bytes()).hexdigest(),'model':'zero-redshift solar MEKAL, wabs attenuation','kT_keV':float(np.exp(best.x[2])),'NH_1e21':float(best.x[1]),'fitted_source_counts':float(np.exp(best.x[0])),'deviance':float(best.fun),'temperature_bounds_keV':[.02,5.0],'optimizer_converged':bool(best.success)}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('spectrum',type=Path);p.add_argument('mekal_mod',type=Path)
    p.add_argument('--fixed-nh',type=float,help='Freeze column in units of 1e21 cm^-2')
    a=p.parse_args();print(json.dumps(run(a.spectrum,a.mekal_mod,a.fixed_nh),indent=2))
