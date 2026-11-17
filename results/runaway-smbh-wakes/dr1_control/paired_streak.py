"""Measure the published DR1 development streak on both visits at the same sky position."""
import argparse
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / 'euclid-foundation-eligibility' / 'followup'
sys.path.insert(0, str(HERE.parent.parent / 'euclid_q1_dark_lenses' / 'arc_search'))
from search import get_range

def load(entry, center, cache):
    key = entry['url'].split('?')[0].split('/ibe/data/euclid/')[1]
    url = 'https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/' + key
    target = cache / (entry['id'] + '-' + entry['kind'] + '.npz')
    receipt = target.with_suffix('.json')
    if not target.exists():
        header_bytes = get_range(url, 0, 28799)
        # Find the first true END card before image data.
        endcard = next(i for i in range(0,len(header_bytes),80) if header_bytes[i:i+8]==b'END     ')
        offset = ((endcard + 80 + 2879)//2880)*2880
        header = fits.Header.fromstring(header_bytes[:offset].decode('ascii'), sep='')
        x, y = WCS(header).world_to_pixel_values(*center)
        xc, yc = round(float(x)), round(float(y))
        w, h = header['NAXIS1'], header['NAXIS2']
        assert min(xc,yc,w-xc,h-yc)>100 and abs(header['BITPIX'])==32
        first, last = offset+4*w*(yc-80), offset+4*w*(yc+81)-1
        data = get_range(url, first,last)
        dtype = '>f4' if header['BITPIX']==-32 else '>i4'
        image = np.frombuffer(data,dtype).reshape(161,w)[:,xc-80:xc+81]
        image = image.astype(float)*header.get('BSCALE',1)+header.get('BZERO',0)
        np.savez_compressed(target,image=image,center=np.array([float(x)-xc+80,float(y)-yc+80]))
        record=dict(url=url,id=entry['id'],kind=entry['kind'],raw_header_sha256=hashlib.sha256(header_bytes[:offset]).hexdigest(),
                    bytes_range=[first,last],raw_range_sha256=hashlib.sha256(data).hexdigest(),x_center=xc,y_center=yc)
        receipt.write_text(json.dumps(record,indent=2)+'\n')
    with np.load(target) as z:
        return entry['id'],entry['kind'],z['image'],z['center'],json.loads(receipt.read_text())

def contrast(image, valid, x, y, angle):
    # Drop the central arcsecond to keep a coincident point source out of the test.
    t=np.r_[np.arange(-44,-6),np.arange(7,45)]
    u=np.r_[-7,-6,-5,-1,0,1,5,6,7]
    xx=x+t[None,:]*np.cos(angle)-u[:,None]*np.sin(angle)
    yy=y+t[None,:]*np.sin(angle)+u[:,None]*np.cos(angle)
    ok=map_coordinates(valid.astype('f4'),[yy,xx],order=1)>0.999
    if 1-ok.mean()>0.10:return np.nan
    a=map_coordinates(image,[yy,xx],order=1)
    a[~ok]=np.nan
    return float(np.nanmedian(np.nanmean(a[3:6],axis=0)-.5*(np.nanmean(a[:3],axis=0)+np.nanmean(a[6:],axis=0))))

def run(source,cache):
    cache.mkdir(parents=True,exist_ok=True)
    prior=json.loads((source/'pilot-v6/trigger-coverage.json').read_text())[0]
    assert prior['x']==2810 and prior['y']==782
    center=(prior['ra_deg'],prior['dec_deg'])
    entries=json.loads((source/'inputs/cutout-requests.json').read_text())
    with ThreadPoolExecutor(4) as p: planes=list(p.map(lambda e:load(e,center,cache), entries[:6]))
    visits={}
    for name,kind,im,c,receipt in planes: visits.setdefault(name,{'center':c})[kind]=im
    angles=np.deg2rad(np.arange(180))
    rows=[]; locked=None
    for name,d in sorted(visits.items()):
        # Keep the Euclid source-detection flag and check RMS coverage separately.
        valid=np.isfinite(d['science']) & (d['rms']>0) & (d['rms']<1e8) & np.isfinite(d['rms'])
        x,y=d['center']
        stats=np.array([contrast(d['science'],valid,x,y,a) for a in angles])
        if locked is None: locked=int(np.nanargmax(stats))
        baselines=[contrast(d['science'],valid,x+dx,y+dy,angles[locked]) for dx in range(-32,33,8) for dy in range(-32,33,8)
                   if abs(-dx*np.sin(angles[locked])+dy*np.cos(angles[locked]))>10]
        b=np.array(baselines); b=b[np.isfinite(b)]; median=np.median(b); mad=1.4826*np.median(abs(b-median))
        rows.append(dict(visit=name,angle_locked_deg=locked,contrast=stats[locked],n_parallel_comparisons=len(b),comparison_median=median,
                         comparison_mad=mad,contrast_in_mad=(stats[locked]-median)/mad))
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.55), layout='constrained')
    for i, ((name,d), ax) in enumerate(zip(sorted(visits.items()), axes)):
        x,y=d['center']
        arr=d['science']
        ax.imshow(np.arcsinh(arr/.0018),origin='lower',cmap='gray',vmin=-1.7,vmax=3.7,
                  extent=[(-x-.5)*.1,(arr.shape[1]-x-.5)*.1,(-y-.5)*.1,(arr.shape[0]-y-.5)*.1])
        for sign in (1,-1):
            t=np.array([4.9,6.35])*sign
            ax.plot(np.cos(angles[locked])*t,np.sin(angles[locked])*t,color='#ffb82c',lw=1.3)
        ax.set_title(f"Visit {i+1} · contrast {rows[i]['contrast_in_mad']:.1f} MAD",fontsize=10)
        ax.set_xlabel('Mosaic x offset (arcsec)')
    axes[0].set_ylabel('Mosaic y offset (arcsec)')
    fig.savefig(HERE/'paired_vis.png',dpi=145)
    plt.close(fig)
    (HERE/'measurement.json').write_text(json.dumps(dict(sky_icrs_deg=center,visits=rows,inputs=[r[-1] for r in planes]),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=SOURCE)
    p.add_argument('--cache',type=Path,required=True)
    a=p.parse_args();run(a.source,a.cache)
