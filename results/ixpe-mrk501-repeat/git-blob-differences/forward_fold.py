"""Exploratory response-folded baseline check, no holdout support.

Power-law shared spectrum and constant q/u across pieces. Uses global exposure
ratio and on-axis responses, so output is a sensitivity check, not a prediction.
"""
import json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy import units as u
from scipy.optimize import least_squares
from development import ROOT,CAL,CENTER,region,particle_keep
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.irf.modf import xModulationFactor

BINS=np.arange(2.,9.,1.)

def responses(du,date):
    with fits.open(CAL/'rmf'/f'ixpe_d{du}_obssim{date}_v013.rmf') as h:
        d=h['MATRIX'].data; true=(d['ENERG_LO']+d['ENERG_HI'])/2
        widths=np.asarray(d['ENERG_HI']-d['ENERG_LO'],float)
        r=np.array(list(d['MATRIX']),float); eb=h['EBOUNDS'].data;rec=(eb['E_MIN']+eb['E_MAX'])/2
        if not (np.all(d['N_GRP']==1) and np.all(d['F_CHAN']==0) and np.all(d['N_CHAN']==375) and np.array_equal(eb['CHANNEL'],np.arange(375))):raise ValueError('unsupported RMF channel groups')
        if np.any(widths<=0):raise ValueError('invalid energy widths')
    if r.shape!=(len(true),375):raise ValueError('unexpected rmf')
    a=xEffectiveArea(str(CAL/'arf'/f'ixpe_d{du}_obssim{date}_v013.arf'))(true)
    m=xModulationFactor(str(CAL/'modfact'/f'ixpe_d{du}_obssim{date}_mfact_v013.fits'))(true)
    transfer=np.column_stack([r[:,(rec>=lo)&(rec<hi)].sum(axis=1) for lo,hi in zip(BINS[:-1],BINS[1:])])
    return true,a,m,transfer,widths

def data(du):
    with fits.open(ROOT/f'auxil/ixpe05004899_det{du}_expmap2_v01.fits.gz') as h:
        w=WCS(h[0].header);a=np.asarray(h[0].data,float);yy,xx=np.indices(a.shape);ra,dec=w.all_pix2world(xx,yy,0);s,b=region(ra,dec);scale=a[s].sum()/a[b].sum()
    with fits.open(ROOT/f'event_l2/ixpe05004899_det{du}_evt2_v01.fits.gz') as h:
        hdr=h['EVENTS'].header; d=h['EVENTS'].data
        t=np.asarray(d['TIME']);q=np.asarray(d['Q']);v=np.asarray(d['U']);pi=np.asarray(d['PI']);e=pi*0.04+0.02
        ew=w.deepcopy();ew.wcs.crpix=[hdr['TCRPX7'],hdr['TCRPX8']];ew.wcs.crval=[hdr['TCRVL7'],hdr['TCRVL8']];ew.wcs.cdelt=[hdr['TCDLT7'],hdr['TCDLT8']]
        ra,dec=ew.all_pix2world(d['X'],d['Y'],1);s,b=region(ra,dec)
        ref=Time(hdr['MJDREFI']+hdr['MJDREFF'],format='mjd',scale='tt');boundary=(Time('2026-07-01T12:00:00',scale='utc').tt-ref).to_value(u.s)
        parent=np.zeros(len(t),np.int8);membership=np.zeros(len(t),np.int8)
        for seg in range(1,6):
            with fits.open(ROOT/f'hk/ixpe0500480{seg}_det{du}_gti_v01.fits.gz') as g:
                a=np.asarray(g['GTI'].data['START']);z=np.asarray(g['GTI'].data['STOP']);k=np.searchsorted(a,t,side='right')-1;mask=(k>=0)&(t<=z[np.maximum(k,0)]);parent[mask]=seg;membership+=mask
        if np.any(membership!=1):raise ValueError('invalid GTI membership')
        chunks=[]
        for seg in range(1,6):
            for date,tm in [('20260101',t<boundary),('20260701',t>=boundary)]:
                pp=(parent==seg)&tm&particle_keep(hdr,d['STATUS'])
                if np.count_nonzero(pp&s)==0:continue
                obs=[]; whiten=[]
                for lo,hi in zip(BINS[:-1],BINS[1:]):
                    bsel=pp&b&(e>=lo)&(e<hi);ssel=pp&s&(e>=lo)&(e<hi)
                    def vec(mask):return np.column_stack((np.ones(mask.sum()),q[mask],v[mask]))
                    sv=vec(ssel);bv=vec(bsel);z=sv.sum(axis=0)-scale*bv.sum(axis=0);c=sv.T@sv+scale**2*(bv.T@bv)
                    try:l=np.linalg.cholesky(c)
                    except np.linalg.LinAlgError:raise ValueError('insufficient covariance')
                    obs.append(z);whiten.append(np.linalg.inv(l))
                chunks.append({'du':du,'segment':seg,'date':date,'obs':np.array(obs),'whiten':np.array(whiten),'response':responses(du,date)})
        return chunks

def fit(chunks):
    def residual(p):
        gamma,q,u=p[:3];out=[]
        for i,c in enumerate(chunks):
            e,area,mu,r,widths=c['response'];flux=e**(-gamma)*area*widths
            pred_i=flux@r;pred_p=(flux*mu)@r
            model=np.column_stack((pred_i,q*pred_p,u*pred_p))*np.exp(p[3+i]);out.extend(np.einsum('bij,bj->bi',c['whiten'],c['obs']-model).ravel())
        return np.array(out)
    p0=np.array([2.5,0.,0.]+[4.]*len(chunks));lo=[0.,-1.,-1.]+[-20.]*len(chunks);hi=[6.,1.,1.]+[20.]*len(chunks)
    fit=least_squares(residual,p0,bounds=(lo,hi),max_nfev=300)
    if not fit.success:raise ValueError(fit.message)
    cov=np.linalg.inv(fit.jac.T@fit.jac)
    return {'gamma':fit.x[0],'q':fit.x[1],'u':fit.x[2],'cov_qu':cov[1:3,1:3].tolist(),'statistic':float(fit.fun@fit.fun),'dof':len(fit.fun)-len(fit.x),'pieces':len(chunks),'normalizations':np.exp(fit.x[3:]).tolist()}

if __name__=='__main__':print(json.dumps(fit(data(1)+data(3)),indent=2))
