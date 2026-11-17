"""Native W2 matched-filter differences at propagated positions with interpolated-grid injections."""
import argparse, csv, hashlib, json, warnings
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy.ndimage import gaussian_filter, maximum_filter, map_coordinates
from scipy.spatial import cKDTree
from acquire import HERE,act_rows,act_sky
from propagate import family,sky,epoch
from other_fields import position
from candidates import circular_family, other_sky

SIGMA = 6.4 / (2.355 * 2.75)
RNG_SEED = 217409

def pix(w,a,d):
    x,y=w.all_world2pix(a,d,0,quiet=True);return np.array([y,x])

def tangent(track,e,n):
    a,d=track;return a+e/(3600*np.cos(np.radians(d))),d+n/3600

def mad(x):
    m=np.nanmedian(x);return m,1.4826*np.nanmedian(np.abs(x-m))

class Frame:
    def __init__(self,p,folder,sigma=SIGMA):
        self.sigma=sigma
        self.visit=p['visit'];self.id=p['id'];self.mjd=p['mjd'];self.t=epoch(Time(p['mjd'],format='mjd',scale='utc').isot)
        with fits.open(folder/(self.id+'_int.fits')) as f:
            h=f[0].header;im=f[0].data.astype('f4');self.w=WCS(h,naxis=2);self.zp=float(h['MAGZP'])
        with fits.open(folder/(self.id+'_msk.fits')) as f:
            self.mw=WCS(f[0].header,naxis=2);mask=(f[0].data.astype('i8') & 1504706175)!=0
        # IBE cutouts and original frames differ only by integer native-pixel offsets.
        cy,cx=np.array(im.shape)//2; xy=np.array([[cx,cy],[cx+100,cy],[cx,cy+100]])
        shift=self.mw.all_world2pix(self.w.all_pix2world(xy,0),0,quiet=True)-xy
        if np.max(np.abs(shift-np.rint(shift[0])))>0.02:raise ValueError('Mask WCS mismatch')
        dx,dy=np.rint(shift[0]).astype(int); yi,xi=np.indices(im.shape,dtype='i2');xi+=dx;yi+=dy
        inside=(xi>=0)&(yi>=0)&(xi<mask.shape[1])&(yi<mask.shape[0]);good=np.zeros(im.shape,dtype=bool)
        good[inside]=~mask[yi[inside],xi[inside]];self.good=good&np.isfinite(im)
        self.zscale=10**(0.4*(20-self.zp));image=np.where(self.good,im,0.)
        self.small=gaussian_filter(self.good.astype('f4'),self.sigma,mode='constant');self.large=gaussian_filter(self.good.astype('f4'),3.,mode='constant')
        self.im=(gaussian_filter(image,self.sigma,mode='constant')/np.maximum(self.small,.001)-gaussian_filter(image,3.,mode='constant')/np.maximum(self.large,.001))*self.zscale
        self.invalid=(self.small<.80)|(self.large<.75)

    def sample(self,track):
        coords=pix(self.w,*track)
        vals=map_coordinates(self.im,coords,order=1,mode='constant',cval=np.nan)
        bad=map_coordinates(self.invalid.astype('f4'),coords,order=1,mode='constant',cval=1)>.01
        vals[bad]=np.nan;return vals

def positions(name,frames,refined=False):
    inputs=json.loads((HERE.parent/'orbits'/'orbit_input.json').read_text()); rng=np.random.default_rng(RNG_SEED)
    if name in ('Phan','Phan-W1'):
        states=family(inputs,*( (28,54,94) if refined else (17,41,73)))
        if refined:
            idx=rng.choice(len(states[0]),size=256,replace=False);states=tuple(x[idx] for x in states[:3])+(states[3],)
        return [sky(states,f.t) for f in frames]
    # Local common positional error only; ACT velocity/distance errors are uncalibrated.
    if refined: e,n=rng.uniform(-59,59,(2,256))
    else:
        e,n=np.meshgrid(np.arange(-60,61,4),np.arange(-60,61,4));e,n=e.ravel(),n.ravel()
    if name.startswith('ACT'):
        return [tangent(act_sky(act_rows()[name],f.t),e,n) for f in frames]
    return [tangent(position(name,f.t),e,n) for f in frames]

def offsets(track,count=512):
    rng=np.random.default_rng(RNG_SEED+1);angle=rng.uniform(0,2*np.pi,count);radius=rng.uniform(210,410,count)
    idx=rng.integers(0,len(track[0][0]),count)
    return [tangent((a[idx],d[idx]),radius*np.cos(angle),radius*np.sin(angle)) for a,d in track]

def gaussian_stamp_response(f,template,truth):
    """Filter unit-total synthetic sources on native pixels and interpolate at trial locations.

    Local signal stamps can be filtered separately because the detector is linear until visit combination.
    Each injected source is a sampled Gaussian normalized to unit integrated counts on that exposure.
    """
    q=pix(f.w,*template);p=pix(f.w,*truth);out=np.zeros(q.shape[1]); rad=7
    yy,xx=np.mgrid[-rad:rad+1,-rad:rad+1];amps=[]
    # Nonzero kernel response only for nearby sky positions.
    for i in np.flatnonzero(np.max(np.abs(q-p),axis=0)<18):
        cy,cx=np.rint(p[:,i]);dy,dx=p[:,i]-[cy,cx]
        stamp=np.exp(-((xx-dx)**2+(yy-dy)**2)/(2*f.sigma**2));stamp/=stamp.sum()
        cyi,cxi=int(cy),int(cx);half=rad+14
        if cyi-half<0 or cxi-half<0 or cyi+half>=f.im.shape[0] or cxi+half>=f.im.shape[1]:continue
        region=np.s_[cyi-half:cyi+half+1,cxi-half:cxi+half+1];pad=np.pad(stamp,14)*f.good[region]
        filt=gaussian_filter(pad,f.sigma)/np.maximum(f.small[region],.001)-gaussian_filter(pad,3.)/np.maximum(f.large[region],.001)
        pos=q[:,i]-[cy,cx]+rad+14;out[i]=map_coordinates(filt,pos[:,None],order=1,mode='constant')[0]
    return out

def main_search(name,cache,out,circular_grid=False):
    rec=json.loads((HERE/'receipts'/(name+'.json')).read_text());folder=cache/name
    band=rec.get('band','W2');sigma=(6.1 if band=='W1' else 6.4)/(2.355*2.75)
    frames=[Frame(p,folder,sigma) for p in rec['products'] if 'int_sha256' in p and 'msk_sha256' in p];visits=sorted(set(f.visit for f in frames))
    min_frames=5 if rec.get('extended_selection') else 2
    visits=[v for v in visits if sum(f.visit==v for f in frames)>=min_frames]
    if len(visits)<2:raise ValueError('fewer than two visits')
    anchor=positions(name,frames);controls=offsets(anchor)
    if circular_grid:
        states=circular_family(name.split('-h')[0]);obs=[other_sky(states,f.t) for f in frames];fake=obs
    else:
        obs=anchor;fake=positions(name,frames,True)
    # Choose disjoint visit for static sky, using predicted displacements alone.
    which={v:np.flatnonzero([f.visit==v for f in frames]) for v in visits};middle={v:which[v][len(which[v])//2] for v in visits}
    def sep(i,j):return np.hypot((np.median(anchor[i][0])-np.median(anchor[j][0]))*np.cos(np.radians(np.median(anchor[i][1]))),np.median(anchor[i][1])-np.median(anchor[j][1]))
    ref={v:max((other for other in visits if other!=v),key=lambda w:sep(middle[v],middle[w])) for v in visits}
    usable=[v for v in visits if len(which[v])>=2 and len(which[ref[v]])>=2 and sep(middle[v],middle[ref[v]])>0.025]
    if len(usable)<2:raise ValueError('fewer than two independent usable visit epochs')
    def features(track):return np.array([z for v in usable for z in (track[middle[v]][0]*np.cos(np.radians(np.median(obs[middle[v]][1]))),track[middle[v]][1])]).T
    dist,near=cKDTree(features(obs)).query(features(fake));best=[(a[near],d[near]) for a,d in obs]
    def raw_eval(tracks,inj=False):
        signal=[]
        for v in usable:
            items=[]
            for i in which[v]:
                track=tracks[i];a=frames[i].sample(track); rr=np.array([frames[j].sample(track) for j in which[ref[v]]]);
                ar=gaussian_stamp_response(frames[i],track,fake[i]) if inj else 0
                br=np.array([gaussian_stamp_response(frames[j],track,fake[j]) for j in which[ref[v]]]) if inj else 0
                items.append((a,rr,ar,br))
            signal.append(items)
        return signal
    def aggregate(raw,flux=0.):
        visitsv=[]
        for items in raw:
            ss=[]
            for a,rr,ar,br in items:
                reference=rr+br*flux
                rv=np.nanmedian(reference,axis=0);rv[np.sum(np.isfinite(reference),axis=0)<2]=np.nan
                ss.append(a+flux*ar-rv)
            s=np.array(ss);x=np.nanmedian(s,axis=0);x[np.isfinite(s).sum(axis=0)<2]=np.nan;visitsv.append(x)
        return np.array(visitsv)
    null=aggregate(raw_eval(controls)); center,scale=np.array([mad(n) for n in null]).T
    if np.any(~np.isfinite(scale)) or np.any(scale<=0):raise ValueError('no control variance')
    def score(x):
        z=(x-center[:,None])/scale[:,None]; eligible=np.all(np.isfinite(z),axis=0); score=np.mean(z,axis=0)
        found=eligible&(np.min(z,axis=0)>=3)&(score>=5); return found,eligible,score,z
    sc,el,val,z=score(aggregate(raw_eval(obs))); nc,ne,_,_=score(null)
    rawfake=raw_eval(best,True);basefake=score(aggregate(rawfake))[0];results=[];indices={}
    for mag in np.arange(12.,19.01,.5):
        detected,eligible,_,_=score(aggregate(rawfake,10**(-.4*(mag-20))));results.append([mag,int((detected&~basefake).sum()),len(near),int(eligible.sum()),int(basefake.sum())])
        if circular_grid:indices[str(mag)]=np.flatnonzero(detected&~basefake).tolist()
    high=np.argsort(np.where(np.isfinite(val),val,-np.inf))[::-1][:10]
    summary={'target':name,'band':band,'frame_count':len(frames),'usable_visits':usable,'visit_minimum_acquired_frames':min_frames,'visit_frames':{str(v):[frames[i].id for i in which[v]] for v in visits},'reference_visits':ref,'visit_median_mjd':{str(v):float(np.median([frames[i].mjd for i in which[v]])) for v in visits},'templates':len(sc),'eligible':int(el.sum()),'threshold_crossings':int(sc.sum()),'control_tracks':len(nc),'eligible_controls':int(ne.sum()),'control_crossings':int(nc.sum()),'control_median':center.tolist(),'control_mad_sigma':scale.tolist(),'injection_nearest_track_rms_as_median_max':(np.array([np.median(dist),np.max(dist)])*3600/np.sqrt(len(usable))).tolist(),'injections_columns':[band+'_Vega_gaussian','recovered_after_blank_veto','tested','eligible','preexisting_triggers'],'injections':results,'highest_mean_scores':[{'template':int(k),'mean':round(float(val[k]),4),'per_visit_z':np.round(z[:,k],3).tolist(),'triggered':bool(sc[k]),'radec_at_first_frame':[round(float(obs[0][0][k]),6),round(float(obs[0][1][k]),6)]} for k in high if el[k]]}
    if circular_grid:summary.update(eligible_track_ids=np.flatnonzero(el).tolist(),trigger_track_ids=np.flatnonzero(sc).tolist(),injection_recovery_ids=indices)
    out.mkdir(parents=True,exist_ok=True);(out/(name+'.json')).write_text(json.dumps(summary,indent=1,allow_nan=False)+'\n');print(name,'eligible/total',el.sum(),len(el),'triggers',sc.sum(),'controls',nc.sum(), 'recovery',results[4:11:2],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('targets',nargs='+');p.add_argument('--circular-grid',action='store_true');args=p.parse_args()
    warnings.filterwarnings('ignore',message='All-NaN slice encountered');warnings.filterwarnings('ignore',category=Warning,module='astropy.wcs')
    for name in args.targets:
        try:main_search(name,args.cache,HERE/('circular_results' if args.circular_grid else 'search_results'),args.circular_grid)
        except (ValueError,FileNotFoundError) as e:print(name,'ineligible:',e,flush=True)
