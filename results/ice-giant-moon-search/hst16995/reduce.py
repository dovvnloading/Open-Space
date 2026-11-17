"""Short-exposure HST detector-plane orbit bank and blind image-level plants."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi, optimize
from scipy.stats import beta
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from geometry import AU, ARCSEC, basis, frequency, tangent_basis

DOMAINS = {'Uranus': [(74500, 84000), (110000, 155000)], 'Neptune': [(68000, 110000), (110000, 190000)]}
RING_BANDS = {'Uranus': [(26000, 72000), (86000, 108000)], 'Neptune': [(39000, 65000)]}
NAMES = {'Uranus': (799, 705), 'Neptune': (899, 808)}
LOCAL_Y, LOCAL_X = np.mgrid[-13:14, -13:14]


def read_csv(name):
    with (ROOT/name).open() as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def peak_center(im, xy):
    x, y = np.rint(xy).astype(int)
    if not (16 <= x < im.shape[1]-16 and 16 <= y < im.shape[0]-16):
        return None
    stamp = im[y-13:y+14,x-13:x+14]
    bg = np.median(stamp[np.hypot(LOCAL_X,LOCAL_Y)>10])
    sm = ndi.gaussian_filter(stamp, 1)
    iy, ix = np.unravel_index(sm[10:17,10:17].argmax(), (7,7))
    ox, oy = ix-3, iy-3
    aperture = np.hypot(LOCAL_X-ox,LOCAL_Y-oy) <= 3.5
    wt = np.maximum(0, stamp-bg)*aperture
    if wt.sum()<=0:
        return None
    return np.array([x+(wt*LOCAL_X).sum()/wt.sum(), y+(wt*LOCAL_Y).sum()/wt.sum()])


def smooth_flux(image, good, sigma):
    weights = ndi.gaussian_filter(good.astype('f4'),sigma)
    num = ndi.gaussian_filter(np.where(good,image,0).astype('f4'),sigma)
    return np.divide(num,weights,where=weights>.01,out=np.zeros_like(num))


def filt(image, good):
    return smooth_flux(image, good, 1.1)-smooth_flux(image, good, 6.)


def prep(cache, planet):
    rows = [r for r in read_csv('selection.csv') if r['target']==planet.upper()]
    ephem = {}
    for r in read_csv('ephemerides.csv'):
        ephem.setdefault(r['rootname'],{})[int(r['body'])]=np.array([float(r[k]) for k in ('ra','dec','r_au','delta_au','light_min')])
    hashes = {r['file']:r['sha256'] for r in read_csv('inputs.csv')}
    body, ctrl = NAMES[planet]
    frames, profiles, reg = [], [], []
    for row in rows:
        name, eph = row['rootname'], ephem[row['rootname']]
        path=cache/(name+'_flt.fits')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=hashes[path.name]:
            raise ValueError('FITS hash mismatch: '+name)
        with fits.open(path) as h:
            s=h['SCI']; image=s.data.copy()
            dq=h['DQ'].data.copy()
            w=WCS(s.header,fobj=h)
            t=(h[0].header['EXPSTART']+h[0].header['EXPEND'])/2
            dt=float(h[0].header['EXPTIME'])
            zp=-2.5*np.log10(s.header['PHOTFLAM']*s.header['PHOTPLAM']**2/2.99792458e18)-48.60
            if h[0].header['FILTER']!='F845M' or s.header['BUNIT']!='ELECTRONS': raise ValueError(name)
        planet_pos = eph[body]
        ra,dec,r,delta,lt = planet_pos
        pix0=w.all_world2pix([planet_pos[:2]],0,quiet=True)[0]
        def relative(sky):
            return w.all_world2pix(np.atleast_2d(sky),0,quiet=True)[0]-pix0
        # One local pixel matrix per detector/readout. Elements map sky radians to pixels.
        jc=np.column_stack([relative([ra+1/3600/np.cos(np.deg2rad(dec)),dec]),relative([ra,dec+1/3600])])*ARCSEC
        _,east,north=tangent_basis(ra,dec)
        u,v,_=basis(planet,t+2400000.5)
        project=jc @ np.array([[east@u,east@v],[north@u,north@v]])/(delta*AU)
        yy,xx=np.indices(image.shape)
        limb=image>np.nanpercentile(image,99.98)*.3
        ly,lx=np.where(limb)
        cx=np.array([np.median(lx),np.median(ly)])
        eq_radius={'Neptune':24764.,'Uranus':25559.}[planet]
        expected_area=np.pi*eq_radius**2*abs(np.linalg.det(jc))/(delta*AU)**2
        if not .6 < limb.sum()/expected_area < 1.4:
            reg.append(dict(rootname=name,visit=row['visit'],limb_x=None,limb_y=None,dx=None,dy=None,
                            control_in_mask=0,control_snr=None,AB_zp=zp,status='excluded_disk_area'))
            continue
        positions={k:cx+relative(vv[:2]) for k,vv in eph.items() if k!=body}
        pred=positions[ctrl]
        disk_radius=np.sqrt(limb.sum()/np.pi)
        isolated=np.linalg.norm(pred-cx)>disk_radius+22
        measured=peak_center(image,pred) if isolated else None
        # Limit the registration adjustment to measurements close to the limb prediction.
        delta_c=np.zeros(2)
        if measured is not None and np.linalg.norm(measured-pred)<3:
            delta_c=measured-pred
            coord=[measured[1]+LOCAL_Y,measured[0]+LOCAL_X]
            stamp=ndi.map_coordinates(image,coord,order=1)
            rr=np.hypot(LOCAL_X,LOCAL_Y)
            stamp-=np.median(stamp[(rr>10.5)&(rr<=13)])
            flux=stamp[rr<=10].sum()
            if flux>300:
                profiles.append(stamp/flux)
        c=cx+delta_c
        positions={k:pos+delta_c for k,pos in positions.items()}
        rr=np.hypot(xx-c[0], yy-c[1])
        # Exclude the disk, its immediate halo and bad pixels before filtering.
        geom=(rr>disk_radius+14)&(xx>=8)&(xx<image.shape[1]-8)&(yy>=8)&(yy<image.shape[0]-8)
        eq=np.linalg.inv(project)@np.stack([(xx-c[0]).ravel(),(yy-c[1]).ravel()])
        r_plane=np.hypot(*eq).reshape(image.shape)
        for lo,hi in RING_BANDS[planet]:
            geom &= ~((r_plane>=lo)&(r_plane<=hi))
        good=geom&(dq==0)&np.isfinite(image)
        # Strong cosmic-ray cores must be inconsistent with even an undersampled optical PSF.
        smooth=ndi.median_filter(image,3); sharp=image-smooth
        cosmic=(sharp>90)&(sharp>4*np.maximum(12.,smooth))
        base_good=good.copy()
        good &= ~cosmic
        control_good=good.copy()
        for k,pos in positions.items():
            x,y=pos
            if not (np.isfinite(pos).all() and -50<x<image.shape[1]+50 and -50<y<image.shape[0]+50):
                continue
            radius=30 if k in [801,701,702,703,704,705] else (10 if k in [808,715] else 7)
            dist=np.hypot(xx-x, yy-y)
            good &= dist>radius
            base_good &= dist>radius
            if k!=ctrl: control_good &= dist>radius
        filtered=filt(image,good)
        # Support near a mask edge is excluded; individual flagged pixels may be interpolated.
        support=geom.copy()
        support &= ndi.distance_transform_edt(geom)>3
        for k,pos in positions.items():
            if np.isfinite(pos).all() and -50<pos[0]<560 and -50<pos[1]<560:
                support &= np.hypot(xx-pos[0],yy-pos[1])>(32 if k in [801,701,702,703,704,705] else 12 if k in [808,715] else 9)
        fixed_support=support.copy()
        support &= ndi.gaussian_filter(good.astype('f4'),1.1)>.76
        noise=np.zeros(image.shape,'f4'); zero=np.zeros_like(noise)
        bins=(rr/12).astype(int)
        for b in np.unique(bins):
            use=(bins==b)&support
            if use.sum()<80: continue
            med=np.median(filtered[use]); mad=1.4826*np.median(abs(filtered[use]-med))
            zero[bins==b]=med; noise[bins==b]=max(0.7,mad)
        support &= noise>0
        # Standardize using the unplanted frame throughout.
        cc = positions[ctrl]
        cinfo=[False,0.]
        if measured is not None:
            coord=[[cc[1]],[cc[0]]]
            in_mask=ndi.map_coordinates(control_good.astype(float),coord,order=1)[0]>.75
            val=ndi.map_coordinates(filt(image,control_good),coord, order=1)[0]
            # Estimate local noise outside the masked moon, using the same annulus.
            ix,iy=np.rint(cc).astype(int); ns=noise[iy,ix]
            if ns>0: cinfo=[bool(in_mask),float((val-zero[iy,ix])/ns)]
        reg.append(dict(rootname=name,visit=row['visit'],limb_x=cx[0],limb_y=cx[1],dx=delta_c[0],dy=delta_c[1],
                        control_in_mask=int(cinfo[0]),control_snr=round(cinfo[1],3),AB_zp=zp,status='used'))
        frames.append(dict(c=c,project=project,mid=t-lt/1440,dt=dt,zp=zp,visit=row['visit'],good=good,
                           image=image,base_good=base_good,fixed_support=fixed_support,score=np.where(support,(filtered-zero)/np.where(noise>0,noise,1),0),
                           zero=zero,noise=noise,support=support,ctrl=positions[ctrl],check=positions[715 if planet=='Uranus' else 808], dist=delta, sun=r))
    if not profiles: raise ValueError('No isolated control PSF')
    # Median registered control stamps, then average azimuthally to suppress nearby structures.
    stamp=np.median(profiles,axis=0); rr=np.hypot(LOCAL_X,LOCAL_Y)
    psf=np.zeros_like(stamp)
    for z in np.unique(np.rint(rr*2)):
        at=np.rint(rr*2)==z;psf[at]=max(0,np.median(stamp[at]))
    psf[rr>10]=0
    psf/=psf.sum()
    return frames,reg,psf,len(profiles)


def orbit(frames, planet, a, phase):
    # Central light time has been removed separately for each exposure.
    mid=np.array([f['mid'] for f in frames]); t=(mid-np.mean(mid))*86400
    theta=phase[:,None]+frequency(planet,a)[:,None]*t
    x,y=[],[]
    for j,f in enumerate(frames):
        u=a*np.cos(theta[:,j]);v=a*np.sin(theta[:,j]); mat=f['project']
        x.append(f['c'][0]+mat[0,0]*u+mat[0,1]*v)
        y.append(f['c'][1]+mat[1,0]*u+mat[1,1]*v)
    return np.array(x).T,np.array(y).T


def grid(planet, domain):
    # Fixed physical radial and azimuthal spacings, recorded even for unsupported tracks.
    parts=[]
    for iz,(lo,hi) in enumerate(domain):
        for radius in np.arange(lo+(100. if iz and lo==domain[iz-1][1] else 0),hi+1,100.):
            n=int(np.ceil(2*np.pi*radius/550))
            parts.append(np.column_stack([np.full(n,radius),np.arange(n)*2*np.pi/n]))
    return np.concatenate(parts)


def samples(frames, planet, bank):
    x,y=orbit(frames,planet,bank[:,0],bank[:,1])
    raw=[]; ok=[]
    for j,f in enumerate(frames):
        xy=[y[:,j],x[:,j]]
        raw.append(ndi.map_coordinates(f['score'],xy,order=1,mode='constant',cval=0))
        ok.append(ndi.map_coordinates(f['support'].astype('f4'),xy,order=1,mode='constant',cval=0)>.995)
    return np.array(raw,dtype='f4').T,np.array(ok).T


def score(raw,ok,visits,dropout=False):
    # Cap single-frame excursions before combining; alternate visit IDs form two disjoint groups.
    raw=np.where(ok,np.clip(raw,-7,7),0)
    groups=[]; sigmas=[]; counts=[]; stable=[]
    for i,v in enumerate(sorted(set(visits))):
        sel=visits==v
        n=ok[:,sel].sum(1)
        z=np.where(n>=2,raw[:,sel].sum(1)/np.sqrt(np.maximum(n,1)),0.)
        groups.append(i%2); sigmas.append(z); counts.append(n)
        worst=(raw[:,sel].sum(1)-np.max(np.where(ok[:,sel],raw[:,sel],-100),axis=1))/np.sqrt(np.maximum(1,n-1))
        stable.append(np.where(n>=2,np.minimum(z/4.5,worst/2),0))
    z=np.array(sigmas); n=np.array(counts); group=np.array(groups)
    supported=(n>=2)
    half=[];eligible=np.ones(len(raw),bool)
    for side in [0,1]:
        active=supported[group==side]; eligible &= active.sum(0)>=2
        half.append(z[group==side].sum(0)/np.sqrt(np.maximum(1,active.sum(0))))
    second=np.partition(z,-2,axis=0)[-2]
    rank=np.minimum(np.minimum(half[0],half[1])/5,second/4.5)
    if dropout:
        rank=np.minimum(rank,np.partition(np.array(stable),-2,axis=0)[-2])
    return np.where(eligible,rank,0),eligible


def plant(frames, planet, psf, a, phase, magnitude,rng):
    x,y=orbit(frames,planet,a,phase)
    maps=[]
    for j,f in enumerate(frames):
        image=f['image'].copy()
        good=f['good'].copy()
        flux=10**(.4*(f['zp']-magnitude))*f['dt']
        # Input magnitude is defined by the 10-pixel-radius count, without infinite-aperture correction.
        for k,(cx,cy) in enumerate(zip(x[:,j],y[:,j])):
            ix,iy=np.rint([cx,cy]).astype(int)
            if not (16<=ix<image.shape[1]-16 and 16<=iy<image.shape[0]-16): continue
            # Integrate over start/middle/end along the physical orbit.
            midphase=phase[k]; delta=frequency(planet,a[k])*f['dt']/2
            stamp=np.zeros_like(psf)
            for sign,w in [(-1,1/6),(0,2/3),(1,1/6)]:
                if sign:
                    xx,yy=orbit(frames,planet,a[k:k+1],np.array([midphase+sign*delta]))
                    px,py=xx[0,j],yy[0,j]
                else: px,py=cx,cy
                stamp+=w*ndi.shift(psf,[py-iy,px-ix],order=1,mode='constant')
            image[iy-13:iy+14,ix-13:ix+14]+=rng.poisson(np.maximum(0,flux*stamp))
            window=np.s_[iy-15:iy+16,ix-15:ix+16]
            patch=image[window]
            smooth=ndi.median_filter(patch,3)
            sharp=patch-smooth
            good[window]=f['base_good'][window] & ~((sharp>90)&(sharp>4*np.maximum(12.,smooth)))
        new=filt(image,good)
        support=f['fixed_support']&(f['noise']>0)&(ndi.gaussian_filter(good.astype('f4'),1.1)>.76)
        maps.append(dict(f,support=support,score=np.where(support,(new-f['zero'])/np.where(f['noise']>0,f['noise'],1),0)))
    return maps,x,y


def control_orbit(frames,planet):
    data=np.array([f['check'] for f in frames])
    def residual(v):
        x,y=orbit(frames,planet,np.array([v[0]]),np.array([v[1]]))
        return (np.column_stack((x[0],y[0]))-data).ravel()
    a={'Neptune':117647.,'Uranus':86005.}[planet]
    ph=min(np.linspace(0,2*np.pi,120),key=lambda p:np.linalg.norm(residual([a,p])))
    fit=optimize.least_squares(residual,[a,ph],x_scale=[1e5,1],loss='soft_l1',f_scale=1)
    return dict(body_id=715 if planet=='Uranus' else 808,a_km=round(fit.x[0],1),median_error_px=round(float(np.median(np.hypot(*fit.fun.reshape(-1,2).T))),3))


def run(args):
    out=args.out;out.mkdir(exist_ok=True,parents=True)
    planet=args.planet;frames,reg,psf,npsf=prep(args.cache,planet)
    write_csv(out/'registration.csv',reg)
    np.savetxt(out/'psf.csv',psf,delimiter=',',fmt='%.8f')
    bank=grid(planet,DOMAINS[planet]);visits=np.array([f['visit'] for f in frames])
    ranks=[];vetted=[]; num_eligible=0
    for start in range(0,len(bank),20000):
        raw,ok=samples(frames,planet,bank[start:start+20000]);r,eligible=score(raw,ok,visits)
        ranks.append(r); vetted.append(score(raw,ok,visits,True)[0]);num_eligible+=int(eligible.sum())
    rank=np.concatenate(ranks);validated=np.concatenate(vetted); hits=np.flatnonzero(rank>=1)
    peaks=[dict(a_km=bank[k,0],phase_rad=bank[k,1],rank=round(rank[k],5)) for k in np.argsort(rank)[-40:][::-1]]
    write_csv(out/'highest_scores.csv',peaks)
    summary=dict(planet=planet,filter='F845M',exposures=len(frames),visits=len(set(visits)),grid_tracks=len(bank),
                 supported_tracks=num_eligible,threshold_crossings=len(hits),top_rank=float(rank.max()),dropout_crossings=int((validated>=1).sum()),top_dropout_rank=float(validated.max()),
                 psf_control_stamps=npsf,control_orbit=control_orbit(frames,planet),median_delta_au=np.median([f['dist'] for f in frames]),
                 median_sun_au=np.median([f['sun'] for f in frames]),reference_emission_mjd=float(np.mean([f['mid'] for f in frames])))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.no_plants: return
    rng=np.random.default_rng(1699501 if planet=='Uranus' else 1699502)
    plants=[]
    # Hold the same randomly drawn off-grid orbits fixed across flux levels.
    for zone,(lo,hi) in enumerate(DOMAINS[planet]):
        radii=rng.uniform(lo+110,hi-110,args.count)
        phases=rng.uniform(0,2*np.pi,args.count)
        trialx,trialy=orbit(frames,planet,radii,phases)
        # Only a positive passing bank track with RMS <= 2 pixels over ALL exposures
        # can count. Evaluating all such trials is equivalent for injection matching.
        neighbors=[]
        for k,(a,p) in enumerate(zip(radii,phases)):
            ix=np.flatnonzero((abs(bank[:,0]-a)<2500)&(abs(np.angle(np.exp(1j*(bank[:,1]-p))))<.08))
            x,y=orbit(frames,planet,bank[ix,0],bank[ix,1])
            d=np.sqrt(np.mean((x-trialx[k])**2+(y-trialy[k])**2,axis=1))
            neighbors.append(ix[d<=2])
        # Batch plants only if their PSF+background stamps are disjoint in every exposure.
        pending=list(range(args.count))
        batches=[]
        while pending:
            batch=[pending.pop(0)]
            for k in pending.copy():
                if all(np.min(np.hypot(trialx[k]-trialx[b],trialy[k]-trialy[b]))>68 for b in batch):
                    batch.append(k);pending.remove(k)
            batches.append(batch)
        for mag in args.mags:
            for batch in batches:
                maps,_,_=plant(frames,planet,psf,radii[batch],phases[batch],mag,rng)
                for k in batch:
                    ix=neighbors[k]
                    if len(ix):
                        raw,ok=samples(maps,planet,bank[ix]);rr,elig=score(raw,ok,visits)
                        va=float(np.max(score(raw,ok,visits,True)[0]));vb=float(np.max(validated[ix]));after=float(np.max(rr));before=float(np.max(rank[ix]));supported=bool(elig.any())
                    else: after=before=va=vb=0.;supported=False
                    plants.append(dict(zone=zone,a_km=round(radii[k],4),phase_rad=round(phases[k],8),AB_mag=mag,
                                  supported=int(supported),baseline_rank=round(before,4),planted_rank=round(after,4),
                                  recovered=int(after>=1 and before<1),dropout_recovered=int(va>=1 and vb<1)))
            write_csv(out/'plants.csv',plants)
    stats=[]
    for zone in range(len(DOMAINS[planet])):
        for mag in args.mags:
            rows=[r for r in plants if r['zone']==zone and r['AB_mag']==mag]
            n=len(rows);k=sum(r['dropout_recovered'] for r in rows)
            stats.append(dict(zone=zone,AB_mag=mag,drawn=n,supported=sum(r['supported'] for r in rows),screen_recovered=sum(r['recovered'] for r in rows),dropout_recovered=k,
                        fraction=round(k/n,4),lower_95=round(float(beta.ppf(.025,k,n-k+1)),4) if k else 0.,
                        upper_95=round(float(beta.ppf(.975,k+1,n-k)),4) if k<n else 1.))
    write_csv(out/'completeness.csv',stats)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('planet',choices=list(DOMAINS));ap.add_argument('cache',type=Path);ap.add_argument('out',type=Path)
    ap.add_argument('--count',type=int,default=32)
    ap.add_argument('--mags',nargs='+',type=float,default=[22.,23.,24.,25.])
    ap.add_argument('--no-plants',action='store_true');run(ap.parse_args())
