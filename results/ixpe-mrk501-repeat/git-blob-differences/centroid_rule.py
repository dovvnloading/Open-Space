"""Fixed development-only pooled XY centroid rule for a future method audit.

No polarization columns are used. No holdout file paths are supported.
"""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from development import ROOT,CENTER,particle_keep

def centroid(offsets):
    center=np.zeros(2)
    for _ in range(5):
        delta=offsets-center
        keep=np.sum(delta**2,axis=1)<60**2
        if np.count_nonzero(keep)<100:raise ValueError('centroid insufficient source events')
        center=offsets[keep].mean(axis=0)
        if np.hypot(*center)>50:raise ValueError('centroid too far from catalog')
    return center

def baseline(segment=None,root=ROOT,obsid="05004899"):
    offsets=[];frame=CENTER.skyoffset_frame()
    for du in (1,3):
        with fits.open(root/f'auxil/ixpe{obsid}_det{du}_expmap2_v01.fits.gz') as h:w=WCS(h[0].header)
        with fits.open(root/f'event_l2/ixpe{obsid}_det{du}_evt2_v01.fits.gz') as h:
            hdr=h['EVENTS'].header;d=h['EVENTS'].data;e=d['PI']*.04+.02
            mask=(e>=2)&(e<8)&particle_keep(hdr,d['STATUS'])
            w.wcs.crpix=[hdr['TCRPX7'],hdr['TCRPX8']];w.wcs.crval=[hdr['TCRVL7'],hdr['TCRVL8']];w.wcs.cdelt=[hdr['TCDLT7'],hdr['TCDLT8']]
            if segment is not None:
                with fits.open(ROOT/f'hk/ixpe0500480{segment}_det{du}_gti_v01.fits.gz') as g:
                    starts=np.asarray(g['GTI'].data['START']);stops=np.asarray(g['GTI'].data['STOP']);t=np.asarray(d['TIME']);idx=np.searchsorted(starts,t,side='right')-1;mask &= (idx>=0)&(t<=stops[np.maximum(idx,0)])
            ra,dec=w.all_pix2world(d['X'][mask],d['Y'][mask],1)
            p=SkyCoord(ra*u.deg,dec*u.deg).transform_to(frame)
            offsets.append(np.column_stack((p.lon.arcsec,p.lat.arcsec)))
    delta=centroid(np.concatenate(offsets));pos=SkyCoord(lon=delta[0]*u.arcsec,lat=delta[1]*u.arcsec,frame=frame).icrs
    return {'offset_arcsec':delta.tolist(),'ra_deg':pos.ra.degree,'dec_deg':pos.dec.degree}

if __name__=='__main__':
    print('all',baseline())
    for seg in range(1,6):print(seg,baseline(seg))
