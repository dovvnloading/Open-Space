"""Exploratory 05004899 only. Not an authorized holdout pipeline.

Uses a single aggregate exposure map per DU; parent-specific support and sky vs
particle background differ, so results are diagnostic until independently checked.
"""
import argparse
import json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from pathlib import Path
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.irf.modf import xModulationFactor
from ixpeobssim.irf.ebounds import channel_to_energy

RA, DEC = 253.46756951679163, 39.76016913205556
CENTER = SkyCoord(RA*u.deg, DEC*u.deg)
ROOT = Path('D:/Open-Space/campaigns/ixpe-mrk501-baseline-bb4ce318/raw')
CAL = Path('D:/Open-Space/software/ixpeobssim-33.2.1/ixpeobssim/caldb/ixpe/gpd/cpf')

def require_region_in_map(w, shape, center, radius_arcsec=270):
    """Require full spherical annulus in the pixel footprint, with margin."""
    cx,cy=w.all_world2pix(center.ra.deg,center.dec.deg,0)
    scales=np.linalg.svd(w.pixel_scale_matrix,compute_uv=False)*3600
    if not np.isfinite([cx,cy]).all() or np.any(~np.isfinite(scales)) or np.min(scales)<=0:raise ValueError('invalid map position or pixel scale')
    radius=radius_arcsec/np.min(scales)+1
    if cx-radius<0 or cy-radius<0 or cx+radius>=shape[1] or cy+radius>=shape[0]:raise ValueError('source/background region clipped by map edge')


def particle_keep(header,status):
    """HEASARC status.NE.bxxxxxxx1xxxxxxxx: zero-based FITS bit 7."""
    if header.get('XPFLGBGD') is not True or status.ndim != 2 or status.shape[1]!=16:
        raise ValueError('missing or malformed particle flag')
    return ~np.asarray(status[:,7],dtype=bool)


def estimate(source, background, scale):
    """Poisson event sums and delta-method ratio covariance, fixed scale."""
    z=source.sum(axis=0)-scale*background.sum(axis=0)
    c=source.T@source+scale**2*(background.T@background)
    if z[0] <= 0: raise ValueError('nonpositive net I')
    i,q,v=z
    jac=np.array([[-q/i**2,1/i,0],[-v/i**2,0,1/i]])
    return z, jac@c@jac.T, c

def region(ra,dec,center=CENTER):
    sep=center.separation(SkyCoord(ra*u.deg,dec*u.deg)).arcsec
    return sep<60, (sep>=120)&(sep<270)

def process(du,center=CENTER,root=ROOT,obsid="05004899",parents=("05004801","05004802","05004803","05004804","05004805")):
    exp=root/f'auxil/ixpe{obsid}_det{du}_expmap2_v01.fits.gz'
    with fits.open(exp) as h:
        a=np.asarray(h[0].data,dtype=float);w=WCS(h[0].header); yy,xx=np.indices(a.shape); ra,dec=w.all_pix2world(xx,yy,0)
        require_region_in_map(w,a.shape,center)
        s,b=region(ra,dec,center); es=a[s].sum();eb=a[b].sum()
        if not np.isfinite(es) or not np.isfinite(eb) or es<=0 or eb<=0:raise ValueError('invalid exposure map')
        if np.count_nonzero(a[s]>0)<0.9*np.count_nonzero(s):raise ValueError('source footprint incomplete')
        scale=es/eb
    p=root/f'event_l2/ixpe{obsid}_det{du}_evt2_v01.fits.gz'
    with fits.open(p) as h:
        hdr=h['EVENTS'].header; d=h['EVENTS'].data
        assert hdr['TIMESYS']=='TT' and hdr['MJDREFI']==57754
        ref=Time(hdr['MJDREFI']+hdr['MJDREFF'],format='mjd',scale='tt')
        boundary=(Time('2026-07-01T12:00:00',scale='utc').tt-ref).to_value(u.s)
        start=(Time('2026-01-01T00:00:00',scale='utc').tt-ref).to_value(u.s)
        end=(Time('2027-01-01T00:00:00',scale='utc').tt-ref).to_value(u.s)
        if np.min(d['TIME']) < start or np.max(d['TIME']) >= end:raise ValueError('unsupported response validity')
        x,y=np.asarray(d['X']),np.asarray(d['Y'])
        event_wcs=w.deepcopy()
        event_wcs.wcs.crpix=[hdr['TCRPX7'],hdr['TCRPX8']]
        event_wcs.wcs.crval=[hdr['TCRVL7'],hdr['TCRVL8']]
        event_wcs.wcs.cdelt=[hdr['TCDLT7'],hdr['TCDLT8']]
        ra,dec=event_wcs.all_pix2world(x,y,1)
        s,b=region(ra,dec,center)
        energy=channel_to_energy(np.asarray(d['PI']));t=np.asarray(d['TIME']); q=np.asarray(d['Q']);v=np.asarray(d['U'])
        keep=(energy>=2)&(energy<8)&particle_keep(hdr,d['STATUS'])
        parent=np.zeros(len(t),dtype=np.int8)
        membership=np.zeros(len(t),dtype=np.int8)
        for segment, parent_id in enumerate(parents,1):
            with fits.open(root/f'hk/ixpe{parent_id}_det{du}_gti_v01.fits.gz') as gh:
                starts=np.asarray(gh['GTI'].data['START']); stops=np.asarray(gh['GTI'].data['STOP'])
                if not np.isfinite(starts).all() or not np.isfinite(stops).all() or np.any(stops<=starts) or np.any(np.diff(starts)<0): raise ValueError('invalid parent GTIs')
                idx=np.searchsorted(starts,t,side='right')-1
                match=(idx>=0)&(t<=stops[np.maximum(idx,0)])
                membership+=match; parent[match]=segment
        if np.any(membership != 1): raise ValueError('parent GTI membership not unique')
        source=[];background=[]; segments=[]
        for date, m in [('20260101',t<boundary),('20260701',t>=boundary)]:
            e=energy[m&keep];qq=q[m&keep];vv=v[m&keep]
            arf=xEffectiveArea(str(CAL/'arf'/f'ixpe_d{du}_obssim{date}_v013.arf'))
            mod=xModulationFactor(str(CAL/'modfact'/f'ixpe_d{du}_obssim{date}_mfact_v013.fits'))
            area=arf(e);mu=mod(e)
            if np.any(area<=0) or np.any(mu<=0):raise ValueError('invalid calibration')
            z=np.stack([1/area,qq/(area*mu),vv/(area*mu)],axis=1)
            source.append(z[s[m&keep]]);background.append(z[b[m&keep]])
            for seg, parent_id in enumerate(parents,1):
                sel=parent[m&keep]==seg
                ss=z[sel&s[m&keep]];bb=z[sel&b[m&keep]]
                if len(ss):
                    zz,cc,_=estimate(ss,bb,scale)
                    segments.append({'parent_id':parent_id,'response_date':date,'source_count':len(ss),'background_count':len(bb),'q':float(zz[1]/zz[0]),'u':float(zz[2]/zz[0]),'cov_qu':cc.tolist()})
    src=np.concatenate(source); bg=np.concatenate(background); vec,cov,cov_iqu=estimate(src,bg,scale)
    alternatives={}
    for label, alpha in [('none',0.0),('geometric',60**2/(270**2-120**2)),('twice_exposure_map',2*scale)]:
        v,c,ci=estimate(src,bg,alpha)
        alternatives[label]={'q':float(v[1]/v[0]),'u':float(v[2]/v[0]),'cov_qu':c.tolist(),'stokes_sum':v.tolist(),'cov_iqu':ci.tolist()}
    return {'du':du,'sufficient':{'source_count':len(src),'background_count':len(bg),'source_sum':src.sum(axis=0).tolist(),'background_sum':bg.sum(axis=0).tolist(),'source_outer':(src.T@src).tolist(),'background_outer':(bg.T@bg).tolist(),'scale':float(scale)},'background_sensitivity':alternatives,'segments':segments,'counts_source':len(src),'counts_background':len(bg),'exposure_map_scale':float(scale),'geometric_scale':60**2/(270**2-120**2),'stokes_sum':vec.tolist(),'cov_iqu':cov_iqu.tolist(),'q':float(vec[1]/vec[0]),'u':float(vec[2]/vec[0]),'cov_qu':cov.tolist(),'boundary_met_tt_seconds':float(boundary)}

def primary(center=CENTER,root=ROOT,obsid="05004899",parents=("05004801","05004802","05004803","05004804","05004805")):
    results=[process(1,center,root,obsid,parents),process(3,center,root,obsid,parents)]
    z=sum((np.array(x['stokes_sum']) for x in results))
    c=sum((np.array(x['cov_iqu']) for x in results))
    i,q,v=z; jac=np.array([[-q/i**2,1/i,0],[-v/i**2,0,1/i]])
    alternatives={}
    for label in ('none','geometric','twice_exposure_map'):
        vectors=[x['background_sensitivity'][label] for x in results]
        zz=sum(np.array(x['stokes_sum']) for x in vectors); cc=sum(np.array(x['cov_iqu']) for x in vectors)
        ii,qq,vv=zz; jj=np.array([[-qq/ii**2,1/ii,0],[-vv/ii**2,0,1/ii]])
        alternatives[label]={'q':qq/ii,'u':vv/ii,'cov_qu':(jj@cc@jj.T).tolist()}
    return {'primary_du':[1,3],'q':q/i,'u':v/i,'cov_qu':(jac@c@jac.T).tolist(),'background_sensitivity':alternatives,'du_sufficient':[x['sufficient'] for x in results]}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--du',type=int,choices=(0,1,2,3),required=True);parser.add_argument('--center',nargs=2,type=float,metavar=('RA','DEC'),default=(RA,DEC));args=parser.parse_args()
    center=SkyCoord(args.center[0]*u.deg,args.center[1]*u.deg)
    if args.du:
        print(json.dumps(process(args.du,center),indent=2))
    else:
        print(json.dumps(primary(center),indent=2))
