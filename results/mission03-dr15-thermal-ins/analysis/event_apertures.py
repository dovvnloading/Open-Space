"""Audit local EPIC imaging events at supplied sky coordinates using PPS files."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

def event_wcs(h):
    w=WCS(naxis=2)
    for attr,suffix in [('crpix','CRPX'),('crval','CRVL'),('cdelt','CDLT'),('ctype','CTYP')]:
        setattr(w.wcs,attr,[h['REF'+ax+suffix] for ax in 'XY'])
    return w

def evaluate(path,ra,dec):
    with fits.open(path) as f:
        h=f['EVENTS'].header
        e=f['EVENTS'].data
        w=event_wcs(h)
        x,y=w.all_world2pix([[ra,dec]],1)[0]
        xx=(e['X'].astype(float)-x)*abs(h['REFXCDLT'])*3600
        yy=(e['Y'].astype(float)-y)*abs(h['REFYCDLT'])*3600
        dist=np.hypot(xx,yy)
        pn=h['INSTRUME']=='EPN'
        good=(e['FLAG']==0)&(e['PATTERN']<=(4 if pn else 12))
        # STDGTI extensions contain per-CCD valid acquisition intervals.
        gti=np.zeros(len(e),dtype=bool)
        for hd in f:
            if hd.name.startswith('STDGTI') and len(hd.data):
                ccd=int(hd.name[6:])
                ci=np.flatnonzero(e['CCDNR']==ccd)
                t=e['TIME'][ci]
                for interval in hd.data:
                    gti[ci] |= ((t>=interval['START'])&(t<=interval['STOP']))
        good &= gti
        bands={'soft':[200,1000],'band1':[200,500],'band2':[500,1000],'medium':[1000,2000],'hard':[2000,10000]}
        b={}
        for label,(lo,hi) in bands.items():
            mask=good&(e['PI']>=lo)&(e['PI']<hi)
            b[label]=[int(np.sum(mask&(dist>=left)&(dist<right))) for left,right in [(0,8),(0,15),(0,20),(25,45)]]
        sel=good&(e['PI']>=200)&(e['PI']<1000)&(dist<15)
        pixels=Counter(zip(e['CCDNR'][sel],e['RAWX'][sel],e['RAWY'][sel]))
        cols=Counter(zip(e['CCDNR'][sel],e['RAWX'][sel]))
        times=np.asarray(e['TIME'][sel],dtype=float)
        summary={'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                 'input_ra':ra,'input_dec':dec,'target':f[0].header.get('OBJECT'),
                 'instrument':h['INSTRUME'],'mode':h['SUBMODE'],'filter':h['FILTER'],
                 'date_start':h['DATE-OBS'],'date_end':h['DATE-END'],
                 'raw_events':len(e),'good_events':int(good.sum()),'aperture_annulus_arcsec':['0-8','0-15','0-20','25-45'],
                 'energy_band_events':b,
                 'strongest_soft_pixel': [[int(i) for i in pix]+[int(v)] for pix,v in pixels.most_common(3)],
                 'strongest_soft_column': [[int(i) for i in col]+[int(v)] for col,v in cols.most_common(3)],
                 'soft_time_span_s':float(np.ptp(times)) if len(times) else None,
                 'soft_pattern':{str(k):v for k,v in sorted(Counter(e['PATTERN'][sel].tolist()).items())},
                 'soft_xy_std_arcsec':[float(np.std(a[sel])) for a in (xx,yy)] if len(times) else None}
        return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('ra',type=float)
    p.add_argument('dec',type=float)
    p.add_argument('directory',type=Path)
    a=p.parse_args()
    files=sorted(a.directory.glob('*[MP]IEVLI*.FTZ'))
    out=[evaluate(f,a.ra,a.dec) for f in files]
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
