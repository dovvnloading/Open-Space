"""Screen complementary native channels in independent visit and orbit coordinates."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys
import warnings
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage
from fetch import ROOT
sys.path.insert(0,str(ROOT.parent))
import ephemeris,geometry,limb_center,prepare


def write(name, rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def geometry_for(paths, times, planet):
    for file in paths:
        with np.load(file) as d:
            if str(d['planet'])!=planet or d['mjd'][0]>times[0] or d['mjd'][-1]<times[-1]:continue
        try:return ephemeris.interpolate(file,times),file
        except ValueError as error:
            if 'gap' not in str(error):raise
    raise ValueError(f'Missing geometry for {planet} {times[0]:.5f}')


def frame(file, paths, planet):
    with fits.open(file,memmap=False) as h:
        header=h[0].header;s=h['SCI'];w=WCS(s.header,naxis=2)
        image=s.data.copy();err=h['ERR'].data.copy();dq=h['DQ'].data.copy()
        times=np.linspace(header['EXPSTART'],header['EXPEND'],9)
        ra0,dec0=s.header['MT_RA'],s.header['MT_DEC']
        visit=header['VISIT_ID'];band=s.header.get('PUPIL',header['PUPIL'])
        if band=='CLEAR' or not band.startswith('F'):band=header['FILTER']
    (_,_,ids,pos),eph=geometry_for(paths,times,planet)
    b=pos[0];pixel_scale=3600*np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))
    radius=b[4,2]/(2*pixel_scale)
    x,y=w.all_world2pix(ra0,dec0,0)
    registration='FITS moving-target WCS'
    try:limb=limb_center.measure(file,radius,True)
    except ValueError:limb=None
    if limb and np.hypot(limb['dx'],limb['dy'])<=4 and abs(limb['radius_px']-radius)<4:
        x,y=limb['x'],limb['y'];registration='measured limb'
    ra,dec=w.all_pix2world(x,y,0)
    M=[]
    for e,n in ((1,0),(0,1)):
        u,v=w.all_world2pix(*prepare.inverse_sky(ra,dec,e/geometry.ARCSEC,n/geometry.ARCSEC),0)
        M.append((u-x,v-y))
    M=np.array(M);cen=np.array([x,y]);ny,nx=image.shape
    good=np.isfinite(image)&np.isfinite(err)&(err>0)&((dq&np.uint32(0xffffffff-4))==0)
    yy,xx=np.indices(image.shape,dtype='f4')
    good&=(xx>5)&(yy>5)&(xx<nx-6)&(yy<ny-6)&((xx-x)**2+(yy-y)**2>(radius+8)**2)
    theta=np.linspace(0,2*np.pi,361)
    los,east,north=geometry.tangent_basis(b[4,0],b[4,1])
    u,v,_=geometry.basis(planet,times[4]+2400000.5)
    m=np.array([[east@u,east@v],[north@u,north@v],[los@u,los@v]])
    # Rasterize the same projected ring envelopes in the actual filter's detector frame.
    rays=np.stack([xx-x,yy-y],axis=-1)@np.linalg.inv(M)
    deprojected=rays@np.linalg.inv(m[:2].T)*(b[4,4]*geometry.AU/geometry.ARCSEC)
    radial=np.linalg.norm(deprojected,axis=-1);pad=4*pixel_scale*b[4,4]*geometry.AU/geometry.ARCSEC/max(.12,abs(np.linalg.det(m[:2])))
    for lo,hi in geometry.RINGS[planet]:good&=(radial<lo-pad)|(radial>hi+pad)
    for k,body in enumerate(ids[1:],1):
        off=geometry.sky_xy(pos[k,:,0],pos[k,:,1],b[:,0],b[:,1],1.,[0,0]);off[:,0]*=-1
        xy=off@M+cen
        rad=max(6.,pos[k,:,2].max()/(2*pixel_scale)+5)
        if body in (701,702,703,704,705,801,808):rad=max(rad,1.5/pixel_scale)
        for a,c in xy:
            if -rad<a<nx+rad and -rad<c<ny+rad:good&=(xx-a)**2+(yy-c)**2>rad**2
    def smooth(width):
        cov=ndimage.gaussian_filter(good.astype('f4'),width)
        return ndimage.gaussian_filter(np.where(good,image,0),width)/np.maximum(cov,.001),cov
    small,cov=smooth(1.1);back,_=smooth(6)
    good&=cov>.82
    residual=small-back
    z=np.zeros_like(residual)
    for top in range(0,ny,48):
        for left in range(0,nx,48):
            sub=np.s_[top:top+48,left:left+48];vals=residual[sub][good[sub]]
            if len(vals)<80:good[sub]=False;continue
            med=np.median(vals);mad=1.4826*np.median(abs(vals-med))
            if not mad>0:good[sub]=False;continue
            z[sub]=(residual[sub]-med)/mad
    t=Time(times,format='mjd',scale='utc').tdb.jd-b[:,4]*geometry.AU/geometry.C/86400
    rec=dict(filename=file.name,planet=planet,filter=band,visit=visit,registration=registration,native_x=round(float(x),3),
             native_y=round(float(y),3),good_pixels=int(good.sum()),geometry=eph.name,mid_mjd=round(float(times[4]),9))
    return dict(t=t,M=M,cen=cen,m=m,delta=b[4,4]*geometry.AU,z=z,good=good.astype('f4'),record=rec)


def sample(f,radius,phases,epoch):
    theta=phases[:,None]+geometry.frequency(f['record']['planet'],radius)*(f['t']-epoch)*86400
    u,v=radius*np.cos(theta),radius*np.sin(theta);m=f['m']
    denom=(f['delta']+m[2,0]*u+m[2,1]*v)/geometry.ARCSEC
    xy=(np.stack([(m[0,0]*u+m[0,1]*v)/denom,(m[1,0]*u+m[1,1]*v)/denom],-1)@f['M'])+f['cen']
    x,y=xy[...,0],xy[...,1]
    ok=ndimage.map_coordinates(f['good'],[y,x],order=1,mode='constant')>.97
    val=ndimage.map_coordinates(f['z'],[y,x],order=1,mode='constant')
    return np.where(ok.all(1),val.mean(1),np.nan)


def search(frames, obj, band, epoch):
    # Multiple chips in the same integration contribute at most one independent sample.
    visits=defaultdict(lambda:defaultdict(list))
    for f in frames:
        r=f['record'];visits[r['visit']][r['filename'].rsplit('_',2)[0]].append(f)
    ph=np.deg2rad(np.arange(0,360,1.));hits=[];out=[];byvisit=defaultdict(lambda:[0,0,0]);nvis=len(visits)
    for rad in range(52000 if obj=='Uranus' else 35000,160001,1000):
        scores=[];bounds=[]
        for vid,exps in sorted(visits.items()):
            per=np.array([np.nanmean([sample(f,rad,ph,epoch) for f in tiles],axis=0) for tiles in exps.values()])
            ok=np.isfinite(per);count=ok.sum(0)
            # Observe all paths in one-exposure visits, but require >=2 independent integrations
            # for discovery evidence. Record how many tracks have any geometric support.
            enough=count>=min(2,len(exps))
            z=np.nansum(np.clip(per,-8,8),axis=0)/np.sqrt(np.maximum(1,count))
            pos=(np.sum(per>=2.5,axis=0)>=min(2,len(exps))) & (z>=6) & enough
            neg=(np.sum(per<=-2.5,axis=0)>=min(2,len(exps))) & (z<=-6) & enough
            byvisit[vid][0]+=int(enough.sum());byvisit[vid][1]+=int(pos.sum());byvisit[vid][2]+=int(neg.sum())
            # A one-exposure visit must itself reach 6; retain visit-level dropout information.
            z=np.where((pos|neg),z,0)
            scores.append(z);bounds.append(enough)
        zs=np.array(scores);bs=np.array(bounds)
        n=(bs.sum(0)>=2).sum();p=(zs>=6).sum(0)>=2;q=(zs<=-6).sum(0)>=2
        for j in np.where(p)[0]:hits.append(dict(planet=obj,filter=band,epoch=epoch,radius_km=rad,phase_deg=j,visit_scores=' '.join(f'{v:.2f}' for v in zs[:,j])))
        out.append((n,p.sum(),q.sum()))
    vout=[dict(planet=obj,filter=band,epoch_emission_jd=epoch,visit=v,unique_exposures=len(visits[v]),tracks_with_support=int(s[0]),positive=int(s[1]),inverted=int(s[2])) for v,s in byvisit.items()]
    return vout,dict(planet=obj,filter=band,epoch_emission_jd=epoch,n_visits=nvis,grid=360*len(out),supported_two_visits=int(sum(r[0] for r in out)),
                     positive_two_visits=int(sum(r[1] for r in out)),inverted_two_visits=int(sum(r[2] for r in out))),hits


def run(files, paths):
    rows=list(csv.DictReader((ROOT.parent/'manifests/products.csv').open()))
    lookup={r['productFilename']:r['planet'] for r in rows}
    grouped=defaultdict(list);records=[]
    for file in sorted(files.glob('*_cal.fits')):
        if 'jw02739' in file.name and lookup[file.name]=='Uranus' and 'long' in file.name and fits.getval(file,'FILTER')=='F300M':continue
        f=frame(file,paths,lookup[file.name]);r=f['record'];records.append(r)
        grouped[(r['planet'],r['filter'])].append(f)
    write('registrations.csv',records);allvis=[];epochs=[];hits=[]
    for (obj,band),frames in sorted(grouped.items()):
        frames.sort(key=lambda f:f['t'][4])
        chunks=np.split(np.arange(len(frames)),1+np.flatnonzero(np.diff([f['t'][4] for f in frames])>4))
        for chunk in chunks:
            group=[frames[i] for i in chunk];epoch=float(np.median([f['t'][4] for f in group]))
            v,q,h=search(group,obj,band,epoch);allvis+=v;epochs.append(q);hits+=h
    write('visits.csv',allvis);write('epochs.csv',epochs)
    with (ROOT/'crossings.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=['planet','filter','epoch','radius_km','phase_deg','visit_scores']);w.writeheader();w.writerows(hits)
    print(len(records),len(allvis),len(epochs),len(hits),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('files',type=Path);p.add_argument('ephem',nargs='+',type=Path)
    a=p.parse_args()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        run(a.files,a.ephem)
