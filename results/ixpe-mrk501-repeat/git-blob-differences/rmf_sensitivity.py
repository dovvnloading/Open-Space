"""Response-only dispersion sensitivity, no observational event access.

Assumes a power-law photon spectrum, constant intrinsic q/u, and uncorrected
on-axis packaged responses. Tests bias of event-energy weighting, not a fit.
"""
import json
import numpy as np
from astropy.io import fits
from development import CAL
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.irf.modf import xModulationFactor


def bias(du,date,index):
    rmf=CAL/'rmf'/f'ixpe_d{du}_obssim{date}_v013.rmf'
    with fits.open(rmf) as h:
        d=h['MATRIX'].data
        true=(np.asarray(d['ENERG_LO'],float)+np.asarray(d['ENERG_HI'],float))/2
        matrix=np.array(list(d['MATRIX']),dtype=float)
        if matrix.shape != (len(true),375):raise ValueError('Unexpected RMF shape')
        channels=h['EBOUNDS'].data
        rec=(np.asarray(channels['E_MIN'],float)+np.asarray(channels['E_MAX'],float))/2
    mask=(rec>=2)&(rec<8);rec=rec[mask];matrix=matrix[:,mask]
    arf=xEffectiveArea(str(CAL/'arf'/f'ixpe_d{du}_obssim{date}_v013.arf'))
    mod=xModulationFactor(str(CAL/'modfact'/f'ixpe_d{du}_obssim{date}_mfact_v013.fits'))
    a_true=arf(true);mu_true=mod(true);a_rec=arf(rec);mu_rec=mod(rec)
    if min(a_rec.min(),mu_rec.min()) <= 0:raise ValueError('Invalid response')
    flux=true**(-index); rate=flux*a_true
    i=np.sum(rate[:,None]*matrix/a_rec[None,:])
    q=np.sum((rate*mu_true)[:,None]*matrix/(a_rec*mu_rec)[None,:])
    return float(q/i)

if __name__=='__main__':
    print(json.dumps([{'du':d,'date':date,'photon_index':g,'estimated_over_true_stokes':bias(d,date,g)} for d in (1,3) for date in ('20260101','20260701') for g in (1.5,2.0,2.5,3.0,3.5)],indent=2))
