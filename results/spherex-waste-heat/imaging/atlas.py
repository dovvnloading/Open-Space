"""Locate the WISE emission on native AllWISE Atlas pixels."""
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
from urllib.request import urlopen
from urllib.parse import urlencode
import numpy as np
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'wise'))
from fetch import position, targets
from measure import arc
BANDS={1:6.1,3:6.5,4:12.0}

def cloud_cutout(url, ra, dec):
    from astropy.io import fits
    from astropy.wcs import WCS
    with urlopen(url, timeout=110) as response:
        native = response.read()
    with fits.open(io.BytesIO(native), do_not_scale_image_data=True) as hdus:
        hdu = hdus[0]
        x, y = WCS(hdu.header).celestial.world_to_pixel_values(ra, dec)
        x0, y0 = int(round(float(x))) - 67, int(round(float(y))) - 67
        data = np.array(hdu.data[y0:y0+135, x0:x0+135])
        if data.shape != (135,135):
            raise ValueError('Requested cutout crosses atlas edge')
        header = hdu.header.copy()
        header['CRPIX1'] -= x0
        header['CRPIX2'] -= y0
        out = io.BytesIO()
        fits.PrimaryHDU(data, header).writeto(out)
    return out.getvalue(), hashlib.sha256(native).hexdigest()


def image(t,band):
    archive=ROOT.parent/'wise'/'archive'
    name=t['label']+'_catwise.csv.gz'
    with gzip.open(archive/name,'rt') as f:
        nearest=min(csv.DictReader(f),key=lambda r:arc(r,t,'catwise'))
    coadd=nearest['source_id'].split('_')[0]+'_ac51'
    r,d=position(t,2010.5)
    url=(f'https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/{coadd[:2]}/{coadd[:4]}/{coadd}/{coadd}-w{band}-int-3.fits?'
         +urlencode(dict(center=f'{r:.7f},{d:.7f}',size='135pix')))
    cache=ROOT/'stamp_cache';cache.mkdir(exist_ok=True)
    file=cache/f'{t["label"]}_W{band}.fits'
    provenance = cache/f'{t["label"]}_W{band}.json'
    if file.exists() and provenance.exists():
        data_source = json.loads(provenance.read_text())
        b = file.read_bytes()
        if data_source['sha256_fits'] != hashlib.sha256(b).hexdigest():
            raise ValueError(f'Cached FITS digest mismatch: {file}')
    else:
        cloud = ('https://nasa-irsa-wise.s3.us-west-2.amazonaws.com/wise/allwise/images/p3am_cdd/'
                 + f'{coadd[:2]}/{coadd[:4]}/{coadd}/{coadd}-w{band}-int-3.fits')
        data_source = dict(url=url)
        try:
            with urlopen(url,timeout=65) as response:
                b = response.read()
            if b.startswith(b'\x1f\x8b'):
                b = gzip.decompress(b)
            if not b.startswith(b'SIMPLE  ='):
                raise ValueError('IBE response was not FITS')
        except Exception as exc:
            b, native_hash = cloud_cutout(cloud, r, d)
            data_source = dict(url=cloud, native_sha256=native_hash, cutout_url=url, cutout_error=str(exc))
        data_source['sha256_fits'] = hashlib.sha256(b).hexdigest()
        data_source['utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        file.write_bytes(b)
        provenance.write_text(json.dumps(data_source,indent=2)+'\n')
    h={}
    for i in range(0,len(b),80):
        card=b[i:i+80].decode('ascii')
        if card[:8].strip()=='END':break
        if card[8:10]=='= ':h[card[:8].strip()]=card[10:].split('/')[0].strip().strip("'").strip()
    assert int(h['BITPIX'])==-32 and int(h['NAXIS'])==2 and h['CTYPE1']=='RA---SIN' and h['CTYPE2']=='DEC--SIN' and float(h['CROTA2'])==0
    a=np.frombuffer(b,dtype='>f4',offset=((i+80+2879)//2880)*2880,count=int(h['NAXIS1'])*int(h['NAXIS2'])).reshape(int(h['NAXIS2']),int(h['NAXIS1'])).astype(float)
    a=a*float(h.get('BSCALE',1))+float(h.get('BZERO',0))
    da,delta,d0=np.radians([r-float(h['CRVAL1']),d,float(h['CRVAL2'])])
    xi=np.cos(delta)*np.sin(da)
    eta=np.sin(delta)*np.cos(d0)-np.cos(delta)*np.sin(d0)*np.cos(da)
    x=np.degrees(xi)/float(h['CDELT1'])+float(h['CRPIX1'])-1
    y=np.degrees(eta)/float(h['CDELT2'])+float(h['CRPIX2'])-1
    scale=float(h['CDELT2'])*3600
    assert 18<x<a.shape[1]-18 and 18<y<a.shape[0]-18
    meta=dict(label=t['label'],band=band,gaia_epoch=2010.5,bytes=len(b),**data_source)
    return a,x,y,scale,meta

def fit(a,x,y,scale,fwhm):
    yy,xx=np.indices(a.shape);dx=xx-x;dy=yy-y
    sigma=fwhm/scale/2.355
    # Fit a ramp simultaneously so extended backgrounds do not bias position as strongly.
    rr=np.hypot(dx,dy)
    good=(rr<max(9,2.4*fwhm/scale)) & np.isfinite(a)
    X,Y=dx[good],dy[good]; z=a[good]
    background=np.median(z)
    local=z[np.hypot(X,Y)<3]
    amp0=max(np.max(local)-background,0.001)
    span=7/scale
    def residual(p):
        amp,cx,cy,off,slx,sly=p
        return amp*np.exp(-((X-cx)**2+(Y-cy)**2)/(2*sigma**2))+off+slx*X+sly*Y-z
    p=least_squares(residual,[amp0,0,0,background,0,0],bounds=([0,-span,-span,-np.inf,-np.inf,-np.inf],[np.inf,span,span,np.inf,np.inf,np.inf])).x
    return p

def analyze(args):
    t,band=args
    a,x,y,s,receipt=image(t,band)
    fit1=fit(a,x,y,s,BANDS[band]);fit2=fit(a,x,y,s,BANDS[band]*1.2)
    amp,cx,cy=fit1[:3]
    sep=s*np.hypot(cx,cy)
    east,north=-s*cx,s*cy
    row=dict(label=t['label'],band=band,east_as=round(east,3),north_as=round(north,3),separation_as=round(sep,3),change_wider_psf_as=round(s*np.hypot(fit2[1]-cx,fit2[2]-cy),3),centroid_at_bound=bool(max(abs(east),abs(north))>6.98))
    yy,xx=np.indices(a.shape)
    recovered=[]
    for ux in (-42,-21,0,21,42):
        for uy in (-42,-21,0,21,42):
            if ux==uy==0: continue
            xx0,yy0=x+ux,y+uy
            if not (15<xx0<a.shape[1]-15 and 15<yy0<a.shape[0]-15):continue
            kernel=np.exp(-((xx-xx0)**2+(yy-yy0)**2)/(2*(BANDS[band]/s/2.355)**2))
            off=fit(a+amp*kernel,xx0,yy0,s,BANDS[band])[1:3]
            recovered.append(s*np.hypot(*off))
    row['local_placements']=len(recovered)
    row['placement_p68_as']=round(np.percentile(recovered,68),3) if recovered else ''
    row['placement_p90_as']=round(np.percentile(recovered,90),3) if recovered else ''
    return row,receipt,(a,x,y,s,cx,cy)

def attempt(args):
    try:
        return analyze(args)
    except Exception as error:
        t, band = args
        return None, dict(label=t['label'], band=band, error=repr(error)), None

def run():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases=[(t,b) for t in targets() for b in BANDS]
    with ThreadPoolExecutor(4) as pool: out=list(pool.map(attempt,cases))
    measured=[r for r,_,_ in out if r is not None]
    if not measured:
        raise ValueError(str([q for _,q,_ in out]))
    with (ROOT/'centroids.csv').open('w') as f:
        w=csv.DictWriter(f,measured[0].keys(),lineterminator='\n');w.writeheader();w.writerows(measured)
    (ROOT/'queries.json').write_text(json.dumps([q for _,q,_ in out],indent=2)+'\n')
    fig,axes=plt.subplots(5,6,figsize=(10,9.7))
    for i,(r,q,values) in enumerate(out):
        ax=axes[(i//3)%5,(i//15)*3+i%3]
        if values is None:
            ax.set_title(f"{q['label']}   W{q['band']}",fontsize=9)
            ax.text(.5,.5,'Unavailable',ha='center',va='center',transform=ax.transAxes)
            ax.set(xticks=[],yticks=[])
            continue
        a,x,y,s,cx,cy=values
        q1,q2=np.nanpercentile(a,[15,99.2]);n=max(q2-q1,1e-6)
        ax.imshow(np.arcsinh((a-q1)/n*4),origin='lower',cmap='gray',vmin=-0.1,vmax=2.6,interpolation='nearest')
        ax.plot(x,y,marker='+',ms=10,color='#37e2ee',mew=1.4)
        ax.plot(x+cx,y+cy,marker='o',fillstyle='none',ms=8,color='#ff9f34',mew=1.1)
        ax.set(xlim=(x-16/s,x+16/s),ylim=(y-16/s,y+16/s),xticks=[],yticks=[])
        ax.set_title(f'{r["label"]}   W{r["band"]}',fontsize=9,pad=2)
    fig.suptitle('AllWISE at the Hephaistos positions',fontsize=16,y=.995)
    fig.text(.5,.005,'32″ panels  ·  North up, east left  ·  Cyan cross: Gaia at 2010.5  ·  Orange circle: fitted WISE centroid',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.017,1,.977),pad=.32)
    fig.savefig(ROOT/'atlas.png',dpi=125);plt.close(fig)

if __name__=='__main__':run()
