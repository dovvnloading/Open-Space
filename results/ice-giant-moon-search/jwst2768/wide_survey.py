"""Search wide moving-target mosaics on a common circular equatorial orbit grid."""
import argparse,csv,json
from pathlib import Path
import sys
import numpy as np
from astropy.time import Time
from scipy.ndimage import map_coordinates
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import AU,ARCSEC,C,basis,frequency,tangent_basis,sky_xy
from prepare import inverse_sky
from search import load_npz,file_hash

EPOCH=2460347.42


def visit_inputs(filename):
    d=load_npz(filename);d['time']=(Time(d['mjd'],format='mjd',scale='utc').tdb.jd-d['path_delta']*AU/C/86400-EPOCH)*86400
    d['los'],d['east'],d['north']=tangent_basis(d['path_ra'],d['path_dec'])
    return d


def projected(d,radius,angles):
    u,v,_=basis('Uranus',EPOCH)
    theta=angles[:,None]+frequency('Uranus',radius)*d['time']
    xyz=radius*(np.cos(theta)[:,:,None]*u+np.sin(theta)[:,:,None]*v)
    denom=d['path_delta']*AU+(xyz*d['los']).sum(-1)
    return np.stack(((xyz*d['east']).sum(-1),(xyz*d['north']).sum(-1)),axis=-1)/denom[:,:,None]*ARCSEC


def sample(d,offsets):
    xy=offsets@d['M']+d['center'];coord=[xy[:,:,1],xy[:,:,0]]
    ok=map_coordinates(d['good'].astype('f4',copy=False),coord,order=1,mode='constant',cval=0,prefilter=False)>.9999
    z=map_coordinates(d['z'],coord,order=1,mode='constant',cval=0,prefilter=False)
    # The five projections are drawn from the same coadd. Their scatter
    # cannot be treated as independent realizations of the background.
    return z.mean(1),ok.all(1)


def rank(z,ok):
    n=ok.sum(1);total=(z*ok).sum(1)/np.sqrt(np.maximum(1,n))
    second=np.partition(np.where(ok,z,-np.inf),-2,axis=1)[:,-2]
    return np.where(n>=2,np.minimum(second/5,total/8),0),total


def run(files,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    ds=[visit_inputs(p) for p in sorted(files)]
    if len(ds)!=3 or len({str(d['visit']) for d in ds})!=3:raise ValueError('three independent visit mosaics required')
    for d in ds:d['good']=d['good'].astype('f4') # one conversion per image, not per orbit radius
    ph=np.deg2rad(np.arange(0,360,.25));radii=np.arange(200000,650001,500)
    outrows=[];candidates=[];idx=0
    with (out/'radii.csv').open('w',newline='') as f,(out/'passing.csv').open('w',newline='') as g:
        w=csv.writer(f);w.writerow(['a_km','trials','supported_2visits','supported_3visits','positive','negative','largest_rank'])
        wr=csv.writer(g);wr.writerow(['a_km','phase_deg','pooled_screen_score','v1_score','v2_score','v3_score','ok1','ok2','ok3','strongest_sidereal_repeat','sidereal_test_eligible'])
        for a in radii:
            angles=[projected(d,a,ph) for d in ds]
            measurements=[sample(d,o) for d,o in zip(ds,angles)]
            z,ok=map(lambda xx:np.stack(xx,1),zip(*measurements))
            pos,total=rank(z,ok);neg,_=rank(-z,ok);passed=np.flatnonzero(pos>=1)
            w.writerow([a,len(ph),int((ok.sum(1)>=2).sum()),int((ok.sum(1)==3).sum()),len(passed),int((neg>=1).sum()),pos.max()])
            if len(passed):
                # Compare fixed-sky hypotheses anchored at the source's predicted
                # middle position in each positive visit. Move the static source
                # through each target-tracked coadd with the actual Uranus motion.
                alt=np.full(len(passed),-np.inf)
                avail=np.zeros(len(passed),bool)
                for i,(d,o) in enumerate(zip(ds,angles)):
                    t=len(d['mjd'])//2
                    ra,dec=inverse_sky(d['path_ra'][t],d['path_dec'][t],o[passed,t,0]/ARCSEC,o[passed,t,1]/ARCSEC)
                    scores=[];bounds=[]
                    for e in ds:
                        sid=sky_xy(ra[:,None],dec[:,None],e['path_ra'],e['path_dec'],1.,(0.,0.));sid[:,:,0]*=-1
                        sc,b=sample(e,sid);scores.append(sc);bounds.append(b)
                    b=np.stack(bounds,1);s=np.stack(scores,1)
                    viable=ok[passed,i]&(z[passed,i]>=5)&b[:,i]&(b.sum(1)>=2)
                    aval=np.partition(np.where(b,s,-np.inf),-2,axis=1)[:,-2]
                    alt=np.maximum(alt,np.where(viable,aval,-np.inf));avail|=viable
                for j,k in enumerate(passed):
                    wr.writerow([a,np.rad2deg(ph[k]),total[k],*z[k],*ok[k].astype(int),alt[j] if avail[j] else '',bool(avail[j])])
            if idx%50==0:f.flush();g.flush();print(a,flush=True)
            idx+=1
    (out/'parameters.json').write_text(json.dumps(dict(input_sha256=[file_hash(p) for p in sorted(files)],epoch_tdb_jd=EPOCH,
        radial_grid_km=[200000,650000,500],phase_step_deg=.25,visits=[str(d['visit']) for d in ds],
        visit_cut=5,pooled_cut=8,require_visit_count=2,coadd_samples=5,independent_exposures_per_mosaic=1,
        trial_family='circular equatorial J2/J4',large_moon_halo_arcsec={'Ariel':5,'Umbriel':5,'Titania':6,'Oberon':6,'Miranda':3,'Puck':1.5}),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('out');p.add_argument('maps',nargs='+')
    a=p.parse_args();run(a.maps,a.out)
