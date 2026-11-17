"""Run the merged PR's peak, injection and aperture routines on calibrated LAWD 66 NIRCam mosaics."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.cluster.hierarchy import fclusterdata

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('repo',type=Path)
p.add_argument('cache',type=Path)
args=p.parse_args()
src=args.repo/'results/jwst-nearby-companions'
sys.path.insert(0,str(src/'miri'))
sys.path.insert(0,str(src/'vetting'))
import search as parent
from fetch_inputs import digest
from archival_images import aperture
warnings.simplefilter('ignore',FITSFixedWarning)
here=Path(__file__).resolve().parent
ra0,dec0=257.8632,-14.7999
cosdec=np.cos(np.deg2rad(dec0))
INNER,OUTER=2.,5.
ANGLES=24
RADII=(2.25,3.25,4.5)
RATIOS=(.0001,.0003,.001,.003,.01,.03,.1)
gen=(here/'inputs.json').read_text()
inputs=json.loads(gen)['inputs']
files={r['filename']:r for r in inputs}
for r in inputs:
    if digest(args.cache/r['filename'])!=r['sha256']:raise ValueError(r['filename'])
gaia=list(csv.DictReader((here/'archives/LAWD66_GaiaDR3.csv').open()))
host=next(x for x in gaia if x['source_id']=='4139531467491239680')
ps1=list(csv.DictReader((here/'archives/PS1_response.csv').open()))
def xy(ra,dec):return np.stack([(ra-ra0)*3600*cosdec,(dec-dec0)*3600],axis=-1)
def read_cat(obs):
    t=Table.read(args.cache/(obs+'_cat.ecsv'))
    return t,xy(t['sky_centroid'].ra.deg,t['sky_centroid'].dec.deg)
def host_expected(year):
    dt=year-2016
    return np.array([float(host['ra'])+dt*float(host['pmra'])/3.6e6/cosdec,float(host['dec'])+dt*float(host['pmdec'])/3.6e6])
def prep(path):
    obs=path.name[:-9]
    cat,catxy=read_cat(obs)
    with fits.open(path) as hd:
        s=hd['SCI'];h=hd[0].header
        assert s.header['BUNIT']=='MJy/sr'
        mjd=s.header['MJD-AVG']
        band=h['PUPIL'] if h['PUPIL']!='CLEAR' else h['FILTER']
        w=WCS(s.header,naxis=2)
        year=Time(mjd,format='mjd').jyear
        guess=host_expected(year)
        dist=np.linalg.norm(catxy-xy(*guess),axis=1)
        close=np.where(dist<.28)[0]
        i=close[np.argmax(cat['aper_total_flux'][close])]
        if cat['aper_total_flux'][i]<.0002:raise ValueError((obs,'host identification'))
        hx,hy=map(float,w.world_to_pixel_values(cat['sky_centroid'][i].ra.deg,cat['sky_centroid'][i].dec.deg))
        x,y=int(round(hx)),int(round(hy))
        stamp=s.data[y-6:y+7,x-6:x+7]
        sy,sx=np.unravel_index(np.nanargmax(stamp),stamp.shape)
        hx,hy=parent.refine(ndi.gaussian_filter(stamp.astype(float),.7),sx,sy)+[x-6,y-6]
        center=w.pixel_to_world_values(hx,hy)
        pos=xy(*center)
        step=.031 if band=='F164N' else .063
        half=int(np.ceil(15/step))
        parent.STEP,parent.HALF=step,half
        grid=parent.sky_grid(*center)
        yy,xx=np.indices((2*half+1,2*half+1))
        px,py=w.world_to_pixel_values(*grid.pixel_to_world_values(xx,yy))
        # Interpolate only a bounded native FITS cutout for efficiency.
        xlo,xhi=int(np.floor(px.min()))-2,int(np.ceil(px.max()))+3
        ylo,yhi=int(np.floor(py.min()))-2,int(np.ceil(py.max()))+3
        if min(xlo,ylo)<0 or xhi>s.data.shape[1] or yhi>s.data.shape[0]:raise ValueError(obs)
        coords=[py-ylo,px-xlo]
        img=ndi.map_coordinates(s.data[ylo:yhi,xlo:xhi].astype(float),coords,order=1)
        wt=ndi.map_coordinates(hd['WHT'].data[ylo:yhi,xlo:xhi].astype(float),coords,order=1)
        valid=np.isfinite(img)&(wt>.3*np.nanmedian(wt[wt>0]))
        good=ndi.minimum_filter(valid,size=31,mode='constant')
        img[~valid]=np.nanmedian(img[valid])
        filt=parent.smooth_response(img)
        rad=np.hypot(xx-half,yy-half)*step
        edges=(0.,1.,2.,3.,4.,5.,7.,10.,15.,22.)
        bins=np.clip(np.digitize(rad,edges)-1,0,len(edges)-2)
        med=np.array([np.median(filt[(bins==j)&good]) for j in range(len(edges)-1)])
        sig=np.array([1.4826*np.median(abs(filt[(bins==j)&good]-m)) for j,m in enumerate(med)])
        noise=sig[bins];score=(filt-med[bins])/noise
        positions,ss=parent.peaks(score,good&(rad>=.2)&(rad<13.))
        rel=(positions-half)*step
        # An independent field star supplies the injected shape, avoiding the
        # microlensed field object near the white dwarf’s core.
        cr=np.linalg.norm(catxy-pos,axis=1)
        nn=cKDTree(catxy).query(catxy,k=2)[0][:,1]
        refmask=(cr>5.5)&(cr<12)&(nn>1.1)&~np.array(cat['is_extended'])&(cat['aper_total_flux']>0)
        choices=np.where(refmask)[0]
        kr=next(j for j in sorted(choices,key=lambda j:-cat['aper_total_flux'][j])
                if good[tuple(np.rint(half+(catxy[j]-pos)[::-1]/step).astype(int))])
        fx,fy=half+(catxy[kr]-pos)/step
        dx,dy=np.mgrid[-10:11,-10:11]
        core=ndi.map_coordinates(img,[fy+dx,fx+dy],order=1)
        rr=np.hypot(dx,dy)
        core-=np.median(core[(rr>=8)&(rr<=10)])
        core*=np.clip((9-rr)/2,0,1)*float(cat['aper_total_flux'][i]/cat['aper_total_flux'][kr])
        return dict(obs=obs,band=band,visit=obs.split('_')[0],mjd=mjd,date=h['DATE-OBS'],cat=cat,catxy=catxy,host_ujy=float(cat['aper_total_flux'][i])*1e6,host_label=int(cat['label'][i]),psf_label=int(cat['label'][kr]),score=score,noise=noise,good=good,psf=core,step=step,half=half,img=img,rad=rad,rel=rel,snr=ss,pos=pos,raw_ra=float(center[0]),raw_dec=float(center[1]),pipeline=h['CAL_VER'],crds=h.get('CRDS_CTX'),host_offset=np.linalg.norm(pos-xy(*guess)))

def register(a,ref):
    va=a['catxy'];vr=ref['catxy']
    mask=(np.linalg.norm(va-a['pos'],axis=1)>7)&(np.linalg.norm(va-a['pos'],axis=1)<50)&~np.array(a['cat']['is_extended'])
    v=va[mask];d,k=cKDTree(vr).query(v)
    delta=vr[k[d<.22]]-v[d<.22]
    for _ in range(3):
        shift=np.median(delta,axis=0)
        sep=np.linalg.norm(delta-shift,axis=1)
        delta=delta[sep<=max(.025,3*np.median(sep))]
    if len(delta)<12:raise ValueError((a['obs'],len(delta)))
    a['shift']=np.median(delta,axis=0)
    a['ntie']=len(delta);a['rms']=np.sqrt(np.mean(np.sum((delta-a['shift'])**2,axis=1)))
    a['hostxy']=a['pos']+a['shift']
    a['tree']=cKDTree(a['rel'])
    a['cat_tree']=cKDTree(a['catxy'])

def writer(name,rows):
    with (here/name).open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)

def simulate(a):
    parent.STEP,parent.HALF=a['step'],a['half']
    rows=[]
    for r in RADII:
        for i in range(ANGLES):
            theta=2*np.pi*i/ANGLES+.037*r
            pt=a['half']+r/a['step']*np.array([np.cos(theta),np.sin(theta)])
            x,y=np.rint(pt).astype(int)
            if not a['good'][y,x]:continue
            states=[]
            for loc in parent.injected_positions(a,pt,(0.,)+RATIOS):
                states.append(bool(np.any(np.linalg.norm(loc-pt,axis=1)<=2)))
            for c,found in zip(RATIOS,states[1:]):
                rows.append(dict(obs=a['obs'],filter=a['band'],radius_arcsec=r,pa_deg=(90-np.rad2deg(theta))%360,ratio=c,baseline_peak=int(states[0]),recovered=int(found and not states[0])))
    return rows

def phot(a,rel):
    rp=3.0
    pt=rel/a['step']+a['half']
    gx,gy=np.rint(pt).astype(int)
    if not a['good'][gy,gx]: return float('nan'),float('nan')
    val=aperture(a['img'],*pt,rp)
    host_ap=aperture(a['img'],a['half'],a['half'],rp)
    calibr=a['host_ujy']/host_ap
    masks=a['good']
    controls=[]
    for th in np.linspace(0,2*np.pi,32,endpoint=False):
        d=np.linalg.norm(rel)*np.array([np.cos(th),np.sin(th)])
        cp=d/a['step']+a['half']
        gx,gy=np.rint(cp).astype(int)
        if np.linalg.norm(d-rel)<.8 or not masks[gy,gx]:continue
        controls.append(aperture(a['img'],*cp,rp)*calibr)
    sky=np.nanmedian(controls)
    v=val*calibr-sky
    err=np.hypot(1.4826*np.nanmedian(abs(np.array(controls)-sky)),.05*v)
    return float(v),float(err)

def models():
    p=args.cache/'bobcat.tar.gz'
    r=json.loads((src/'miri/calibration_inputs.json').read_text())[0]
    if not p.exists():
        from urllib.request import urlopen
        p.write_bytes(urlopen(r['url'],timeout=120).read())
    if digest(p)!=r['sha256']:raise ValueError(p)
    rows=[]
    with tarfile.open(p) as tar:
        for mh in ('-0.5','+0.0','+0.5'):
            lines=tar.extractfile('photometry_tables/flux_table_JWST'+mh).read().decode().splitlines()
            j=next(i for i,s in enumerate(lines) if s.strip().startswith('Teff'))
            names=lines[j].split();names=names[names.index('F070W'):]
            grid=np.loadtxt(io.StringIO('\n'.join(s for s in lines[j+1:] if s.strip() and s.strip()[0].isdigit() and '*' not in s)))
            for v in grid:
                f=1000*(10*float(host['parallax'])/1000)**2*10**v[[6+names.index(k) for k in ('F164N','F323N','F277W')]]
                rows.append([float(mh),v[0],v[1],*f])
    return np.array(rows)

allobs=sorted([prep(args.cache/r['filename']) for r in inputs if r['filename'].endswith('i2d.fits')],key=lambda a:(a['mjd'],a['band']))
ref=next(x for x in allobs if x['band']=='F164N')
for a in allobs:register(a,ref)
grid=models()
trials=[];scores=[];matches=[];meta=[]
for a in allobs:
    trials+=simulate(a)
    other=next(x for x in allobs if x['visit']==a['visit'] and x['obs']!=a['obs'])
    radii=np.linalg.norm(a['rel'],axis=1)
    area_search=a['good']&(a['rad']>=INNER)&(a['rad']<=OUTER)
    area_control=a['good']&(a['rad']>=7)&(a['rad']<13)
    ncontrol=np.count_nonzero((radii>=7)&(radii<13))
    meta.append(dict(obs=a['obs'],filter=a['band'],date=a['date'],mjd=a['mjd'],ra=a['raw_ra'],dec=a['raw_dec'],host_catalog_label=a['host_label'],injection_psf_catalog_label=a['psf_label'],host_flux_uJy=a['host_ujy'],host_vs_gaia_pm_only_arcsec=a['host_offset'],east_tie_shift=a['shift'][0],north_tie_shift=a['shift'][1],ties=a['ntie'],tie_rms_arcsec=a['rms'],search_area_arcsec2=area_search.sum()*a['step']**2,search_fraction=area_search.sum()/np.count_nonzero((a['rad']>=INNER)&(a['rad']<=OUTER)),control_area_arcsec2=area_control.sum()*a['step']**2,control_peaks=ncontrol,pipeline=a['pipeline'],crds=a['crds']))
    gaia_ref=xy(np.array([float(g['ra']) for g in gaia]),np.array([float(g['dec']) for g in gaia]))
    gaia_pm=np.array([[float(g['pmra'] or 0),float(g['pmdec'] or 0)] for g in gaia])
    gpos=gaia_ref+gaia_pm*(Time(a['mjd'],format='mjd').jyear-2016)/1000
    ppos=xy(np.array([float(s['raMean']) for s in ps1]),np.array([float(s['decMean']) for s in ps1]))
    host_ps1=np.array([xy(*host_expected(y)) for y in Time([float(s['epochMean']) for s in ps1],format='mjd').jyear])
    for i in np.flatnonzero((radii>=INNER)&(radii<=OUTER)):
        rel=a['rel'][i];sky=a['hostxy']+rel;srcid=a['obs']+f'-p{i:04d}'
        visits=[]
        for b in allobs:
            if b['visit']==a['visit']:continue
            dc,ic=b['tree'].query(rel)
            db,ib=b['tree'].query(sky-b['hostxy'])
            motion=np.linalg.norm(a['hostxy']-b['hostxy'])
            rec=dict(id=srcid,comparison_obs=b['obs'],host_shift_arcsec=motion,fixed_resid_arcsec=db,host_resid_arcsec=dc,fixed_nearest_peak=int(ib),moving_nearest_peak=int(ic))
            matches.append(rec)
            if b['band']==a['band'] and motion>=.22:visits.append((dc,db))
        nstat=sum(db<=.10 and db+.08<dc for dc,db in visits)
        ncom=sum(dc<=.10 and dc+.08<db for dc,db in visits)
        ambiguous=sum(min(dc,db)<=.10 for dc,db in visits)-nstat-ncom
        v1,e1=phot(a,rel)
        v2,e2=phot(other,sky-other['hostxy'])
        iband={'F164N':3,'F323N':4,'F277W':5}
        chi=((grid[:,iband[a['band']]]-v1)/e1)**2+((grid[:,iband[other['band']]]-v2)/e2)**2
        j=np.argmin(chi)
        d_cat,ci=a['cat_tree'].query(sky-a['shift'])
        ds=np.linalg.norm(gpos-sky,axis=1);g=gaia[int(np.argmin(ds))]
        ds2=np.linalg.norm(ppos-sky,axis=1);ip=int(np.argmin(ds2));pp=ps1[ip]
        # Compare the host-relative hypothesis at every PS1 object's own reference epoch.
        dp_move=np.linalg.norm(ppos-host_ps1-rel,axis=1)
        dg_move=np.linalg.norm(gaia_ref-xy(*host_expected(2016.))-rel,axis=1)
        radc=np.linalg.norm(a['rel'],axis=1)
        control_count=np.count_nonzero((radc>=7)&(radc<13)&(a['snr']>=a['snr'][i]))
        scores.append(dict(id=srcid,obs=a['obs'],band=a['band'],east_arcsec=rel[0],north_arcsec=rel[1],ra_registered=ra0+sky[0]/3600/cosdec,dec_registered=dec0+sky[1]/3600,r_arcsec=radii[i],score=a['snr'][i],independent_same_band_tests=len(visits),fixed_wins=nstat,host_wins=ncom,ambiguous=ambiguous,cat_offset_arcsec=d_cat,cat_label=int(a['cat']['label'][ci]),cat_extended=int(a['cat']['is_extended'][ci]),flux_uJy=v1,error_uJy=e1,other_band=other['band'],other_flux_uJy=v2,other_error_uJy=e2,model_chi2=chi[j],model_Teff=grid[j,1] if np.isfinite(chi[j]) else '',model_logg=grid[j,2] if np.isfinite(chi[j]) else '',model_metallicity=grid[j,0] if np.isfinite(chi[j]) else '',control_at_least_as_strong=control_count,expected_same_score_in_annulus=control_count*area_search.sum()/area_control.sum(),gaia_id=g['source_id'],gaia_fixed_arcsec=min(ds),ps1_id=pp['objID'],ps1_fixed_arcsec=min(ds2),ps1_moving_id=ps1[int(np.argmin(dp_move))]['objID'],ps1_moving_arcsec=min(dp_move),gaia_moving_id=gaia[int(np.argmin(dg_move))]['source_id'],gaia_moving_arcsec=min(dg_move)))
    print(a['obs'],len(scores),flush=True)
writer('epochs.csv',meta);writer('peaks.csv',scores);writer('tracks.csv',matches);writer('injections.csv',trials)
curves=[]
for a in allobs:
    for r in RADII:
        for ratio in RATIOS:
            sub=[t for t in trials if (t['obs'],t['radius_arcsec'],t['ratio'])==(a['obs'],r,ratio)]
            curves.append(dict(obs=a['obs'],radius_arcsec=r,ratio=ratio,trials=len(sub),recovered=sum(t['recovered'] for t in sub)))
writer('contrast.csv',curves)

# Test the same fixed/host competition on inserted signals at the widest
# same-passband baseline. Preserve real late-epoch peaks outside the local patch.
joint=[]
for band in ('F164N','F277W','F323N'):
    group=[a for a in allobs if a['band']==band]
    a,b=group[0],group[-1]
    step=a['step'];half=a['half']
    parent.STEP,parent.HALF=step,half
    parent.MATCH,parent.MARGIN=.10,.08
    bgshift=(a['hostxy']-b['hostxy'])/step
    for rad in RADII:
        for i in range(ANGLES):
            theta=2*np.pi*i/ANGLES+.037*rad
            pt=half+rad/step*np.array([np.cos(theta),np.sin(theta)])
            x,y=np.rint(pt).astype(int)
            if not a['good'][y,x] or not b['good'][y,x]:continue
            full=b['rel']/step+half
            # peaks() refines pixel-centered local maxima by up to 0.8 pixel.
            far=full[np.any(abs(full-[x,y])>17.8,axis=1)]
            baseline=False
            for ratio,al,bl in zip((0.,)+RATIOS,parent.injected_positions(a,pt,(0.,)+RATIOS),parent.injected_positions(b,pt,(0.,)+RATIOS)):
                late=np.vstack([far,bl])
                ok=False
                for p0 in al[np.linalg.norm(al-pt,axis=1)<=2.]:
                    state,_,_,_,ic=parent.classify(p0,late,bgshift)
                    if state=='common-motion trigger' and np.linalg.norm(late[ic]-pt)<=2.:
                        ok=True
                if ratio==0:baseline=ok
                else:joint.append(dict(band=band,early=a['obs'],late=b['obs'],radius_arcsec=rad,pa_deg=(90-np.rad2deg(theta))%360,ratio=ratio,baseline_match=int(baseline),recovered=int(ok and not baseline)))
writer('joint_injections.csv',joint)

# Follow the marginal single-epoch peak at both predicted loci, including
# below-threshold photometry and all independently produced pipeline catalogs.
forced=[]
for r in scores:
    if r['fixed_wins'] or r['host_wins']:continue
    a=next(v for v in allobs if v['obs']==r['obs'])
    rel=np.array([r['east_arcsec'],r['north_arcsec']])
    for b in allobs:
        for hypothesis,pred in [('host',rel),('fixed',a['hostxy']+rel-b['hostxy'])]:
            flux,err=phot(b,pred)
            dist,idx=b['cat_tree'].query(b['pos']+pred)
            x,y=np.rint(b['half']+pred/b['step']).astype(int)
            forced.append(dict(id=r['id'],prediction=hypothesis,obs=b['obs'],date=b['date'],band=b['band'],flux_uJy=flux,err_uJy=err,footprint_valid=int(b['good'][y,x]),catalog_nearest_arcsec=dist,catalog_nearest_label=int(b['cat']['label'][idx]),catalog_nearest_flux_uJy=float(b['cat']['aper_total_flux'][idx])*1e6))
if forced:writer('unmatched_forced.csv',forced)

# Consolidate repeat peaks within the search annulus using a 0.10 arcsec
# maximum pairwise distance. Track tests can use independent comparison peaks
# outside that annulus, so a group with one in-annulus row need not be unmatched.
sky_all=xy(np.array([r['ra_registered'] for r in scores]),np.array([r['dec_registered'] for r in scores]))
groups=fclusterdata(sky_all,.10,criterion='distance',method='complete')
sources=[]
for group in sorted(set(groups)):
    rows=[r for r,g in zip(scores,groups) if g==group]
    sources.append(dict(group=int(group),ra_registered_deg=np.mean([r['ra_registered'] for r in rows]),dec_registered_deg=np.mean([r['dec_registered'] for r in rows]),in_annulus_triggers=len(rows),bands=' '.join(sorted({r['band'] for r in rows})),unique_visits=len({r['obs'].split('_')[0] for r in rows}),peak_ids=' '.join(r['id'] for r in rows),fixed_wins=max(r['fixed_wins'] for r in rows),host_wins=max(r['host_wins'] for r in rows),min_gaia_fixed_arcsec=min(r['gaia_fixed_arcsec'] for r in rows),min_ps1_fixed_arcsec=min(r['ps1_fixed_arcsec'] for r in rows)))
writer('sources.csv',sources)
jcurves=[]
for band in ('F164N','F277W','F323N'):
    for radius in RADII:
        for ratio in RATIOS:
            rows=[r for r in joint if (r['band'],r['radius_arcsec'],r['ratio'])==(band,radius,ratio)]
            jcurves.append(dict(band=band,radius_arcsec=radius,ratio=ratio,recovered=sum(r['recovered'] for r in rows),trials=len(rows)))
writer('joint_contrast.csv',jcurves)
