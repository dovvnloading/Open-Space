"""Rebuild XMM compact photons, CCD good times and reference translations from HEASARC and Chandra level-2 events.
Usage: python rebuild_inputs.py --raw-cache CACHE --chandra-evt2 PATH --output DATA --manifest xmm-products.csv
The cache may be empty. Chandra ObsID 13813 can be downloaded from the URL in README.md.
"""
import argparse, hashlib, os, shutil, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter,maximum_filter
from scipy.spatial import cKDTree
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
R=15*(13+29/60+39.94/3600);D=47+12/60+36.6/3600;R2=15*(13+29/60+39.44/3600);D2=47+12/60+43.3/3600
XY2=np.array([(R2-R)*np.cos(np.radians(D))*3600,(D2-D)*3600])
def world(hdu,x,y):
 cols=[hdu.columns.names.index(n)+1 for n in ['X','Y']] if 'X' in hdu.columns.names else [hdu.columns.names.index(n)+1 for n in ['x','y']]
 w=WCS(naxis=2)
 for key,attr in [('TCRVL','crval'),('TCRPX','crpix'),('TCDLT','cdelt'),('TCTYP','ctype')]:setattr(w.wcs,attr,[hdu.header[f'{key}{i}'] for i in cols])
 a,d=w.all_pix2world(x,y,1);return (a-R)*np.cos(np.radians(D))*3600,(d-D)*3600
def chandra():
 f='chandra_ref.csv'
 if Path(f).exists():return pd.read_csv(f)
 with fits.open(CHANDRA) as h:
  e=h['EVENTS'];g=e.data;mask=(g['energy']>700)&(g['energy']<7000)
  x,y=world(e,g['x'][mask],g['y'][mask]);xy=np.stack([x,y],axis=1);tree=cKDTree(xy)
 b=np.arange(-620.,621.,1.)
 H,xr,yr=np.histogram2d(x,y,bins=(b,b));H=gaussian_filter(H,.9);xx,yy=np.nonzero((H==maximum_filter(H,13))&(H>4));refs=[]
 for a,d in zip(xx,yy):
  xy0=np.array([xr[a]+.5,yr[d]+.5]); ix=tree.query_ball_point(xy0,2.5)
  if len(ix)<40:continue
  xy1=xy[ix].mean(axis=0);refs.append((*xy1,len(ix)))
 df=pd.DataFrame(refs,columns=['x','y','n']);df.to_csv(f,index=False);return df
def extract(fname):
 out=Path(fname.name[:14]+'.npz')
 if out.exists():return out
 with fits.open(fname) as h:
  e=h['EVENTS'];hdr=e.header;g=e.data;en=g['PI'].copy();mask=(en>=500)&(en<=7000)&(g['FLAG']==0)&(g['PATTERN']<=(4 if 'PN' in fname.name else 12))
  x,y=world(e,g['X'][mask],g['Y'][mask]);t=hdr['MJDREF']+(g['TIME'][mask]+hdr.get('TIMEZERO',0))/86400;en=en[mask]
  a,d=np.radians([R,D]);v=np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
  times=Time([t.min(),t.max()],format='mjd',scale='tt');delta=(v@get_body_barycentric('earth',times.tdb).xyz.to_value('au'))*499.00478/86400
  t=t+np.interp(t,[t.min(),t.max()],delta)
  keep=(abs(x)<620)&(abs(y)<620)
  np.savez_compressed(out,x=x[keep],y=y[keep],e=en[keep],t=t[keep])
 return out
def register(data,ref,low=900):
 keep=(data['e']>=low);pts=np.column_stack([data['x'][keep],data['y'][keep]]);tree=cKDTree(pts)
 refs=ref.loc[(ref.n>=45)&(np.hypot(ref.x,ref.y)>30)&(np.hypot(ref.x,ref.y)<310)].copy()
 srctree=cKDTree(ref[['x','y']]);regs=[]
 for _,r in refs.iterrows():
  xy=np.array([r.x,r.y]);near=ref.iloc[srctree.query_ball_point(xy,20)]
  if (near.n>r.n*.20).sum()>1:continue
  ix=tree.query_ball_point(xy,17); outer=tree.query_ball_point(xy,25);expect=(len(outer)-len(ix))*(17**2/(25**2-17**2))
  if (len(ix)-expect)<45 or expect>.72*len(ix):continue
  pos=xy.copy()
  for i in range(16):
   sub=pts[tree.query_ball_point(pos,14)]
   weights=np.exp(-(np.linalg.norm(sub-pos,axis=1)/6.8)**2/2)
   pos=np.average(sub,axis=0,weights=weights)
  dx=pos-xy
  if np.linalg.norm(dx)<7:regs.append((*dx,*xy,len(ix)))
 return np.array(regs)

def prepare(fn):
 name=fn.name[:14]; sub=DEST/(name+'-photons.csv.gz')
 if sub.exists(): return
 raw=next(Path('.').glob(name+'*EVLI*FTZ')); e=np.load(fn);x,y,en,t=[e[k] for k in ['x','y','e','t']]; radius=np.hypot(x,y)
 near=(radius<38)&(en>=500)
 if near.sum()<100:return
 reg=register(e,chandra())
 if len(reg)<2:return
 # 500-s bins in source-free wing (galactic point sources contribute a stable pedestal).
 bins=np.arange(t.min(),t.max()+500/86400,500/86400);counts,_=np.histogram(t[(radius>120)&(radius<400)&(en>2500)],bins=bins)
 # Use relative 25th percentile. Fixed rule shared across detectors.
 select=(counts>0)&(counts<=1.8*np.percentile(counts[counts>0],25))
 with fits.open(raw) as h:
  ev=h['EVENTS']; g=ev.data; hdr=ev.header
  stride=max(1,len(g)//500000); sx,sy=world(ev,g['X'][::stride],g['Y'][::stride]); chips=g['CCDNR'][::stride][np.hypot(sx,sy)<12]
  if not len(chips):return
  chip=int(np.bincount(chips).argmax());gti=h[f'STDGTI{chip:02d}'].data
  start=hdr['MJDREF']+(gti['START']+hdr.get('TIMEZERO',0))/86400
  stop=hdr['MJDREF']+(gti['STOP']+hdr.get('TIMEZERO',0))/86400
  a,d=np.radians([R,D]);v=np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
  knots=np.array([start.min(),stop.max()]);time=Time(knots,format='mjd',scale='tt');delta=(v@get_body_barycentric('earth',time.tdb).xyz.to_value('au'))*499.00478/86400
  start=start+np.interp(start,knots,delta);stop=stop+np.interp(stop,knots,delta)
 good=[]
 for a,b in zip(bins[:-1][select],bins[1:][select]):
  lo=np.maximum(start,a); hi=np.minimum(stop,b);good.extend(zip(lo[hi>lo],hi[hi>lo]))
 if not good:return
 good=np.array(good); good=good[good[:,0].argsort()]
 # Overlapping or touching raw GTIs have identical exposure treatment.
 intervals=[]
 for a,b in good:
  if intervals and a<=intervals[-1][1]+1e-10:intervals[-1][1]=max(b,intervals[-1][1])
  else:intervals.append([a,b])
 good=np.array(intervals)
 j=np.clip(np.searchsorted(good[:,0],t[near])-1,0,len(good)-1)
 inside=(t[near]>=good[j,0])&(t[near]<=good[j,1])
 shift=np.median(reg[:,:2],axis=0)
 pd.DataFrame({'x_arcsec':np.round(x[near][inside]-shift[0],4),'y_arcsec':np.round(y[near][inside]-shift[1],4),'energy_ev':en[near][inside],'mjd_tdb_geo':np.round(t[near][inside],8)}).to_csv(sub,index=False)
 pd.DataFrame(good,columns=['start_mjd','stop_mjd']).to_csv(DEST/(name+'-gti.csv'),index=False)
 pd.DataFrame(reg,columns=['dx','dy','ref_x','ref_y','events']).to_csv(DEST/(name+'-references.csv'),index=False)
 print(name,len(reg),int(inside.sum()),round(np.sum(good[:,1]-good[:,0])*86.4,1),flush=True)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-cache',type=Path,required=True)
    parser.add_argument('--chandra-evt2',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,default=Path(__file__).with_name('xmm-products.csv'))
    args=parser.parse_args()
    CHANDRA=args.chandra_evt2.resolve()
    DEST=args.output.resolve()
    manifest=pd.read_csv(args.manifest)
    DEST.mkdir(parents=True,exist_ok=True)
    args.raw_cache.mkdir(parents=True,exist_ok=True)
    os.chdir(args.raw_cache)
    for item in manifest.itertuples(index=False):
        raw=Path(item.filename)
        if not raw.exists():
            temp=raw.with_suffix('.part')
            with urllib.request.urlopen(item.url,timeout=240) as response, temp.open('wb') as fout:
                shutil.copyfileobj(response,fout)
            temp.replace(raw)
        with raw.open('rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=item.sha256:
                raise RuntimeError(f'Archive checksum mismatch for {raw}')
        prepare(extract(raw))
