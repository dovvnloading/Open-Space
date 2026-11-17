"""Project the two October 7 hypotheses onto disjoint October 6 F162M data."""
import argparse,csv,json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from link_visits import make_cube,search,trials,evaluate
from independent_vet import save
ROOT=Path(__file__).resolve().parent


def run(args):
    raw=json.loads((ROOT/'independent_vet.json').read_text())['F150W2']
    displaced=json.loads((ROOT/'integration_check/displaced_circular.json').read_text())
    pars=[('screen',raw['radius_km'],raw['phase_deg']),('displaced',displaced['radius_km'],displaced['phase_deg'])]
    audits=[]
    for path in [args.long,args.short]:
        c=make_cube(path,args.ee,False)
        kernels=search.validate_cube(c);cuts=search.Cuts(psf_support=.8);maps=search.prepare(c,kernels,cuts)
        with TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            for name,a,phi in pars:
                # Nearby tracks share the same radius. Keep raw per-frame support.
                phases=phi+np.r_[0,np.arange(-30,0,1),np.arange(1,31,1)]
                el=np.zeros((len(phases),6));el[:,0]=a;el[:,-1]=phases%360
                b=trials(tmp/'t.csv',path,tmp/'b.npz',el,raw['emission_epoch_tdb_jd'])
                n,i=evaluate(c,b,kernels,maps,cuts)
                z=search.snr(n,i);fl=np.divide(n,i,out=np.zeros_like(n),where=i>0)
                for j,obs in enumerate(c['exposure_id']):
                    same=(i[:,j]>0);same[0]=False;far=same&(abs(phases-phi)>=5)
                    nearest=np.linalg.norm(c['known_path'][:,j,:,None,:]-np.stack([b['x_path'][0,j],b['y_path'][0,j]],axis=1)[None,None,:,:],axis=-1).min(axis=(1,2))
                    audits.append(dict(hypothesis=name,exposure_id=obs,visit=c['visit'][j],a_km=round(a,3),phase_deg=round(phi,6),
                                     supported=int(i[0,j]>0),score=round(z[0,j],4),flux_uJy=round(fl[0,j]*1e6,5) if i[0,j]>0 else '',
                                     statistical_error_uJy=round(1e6/np.sqrt(i[0,j]),5) if i[0,j]>0 else '',
                                     local_reference_median_uJy=round(np.median(fl[far,j])*1e6,4) if far.any() else '',
                                     closest_moon=int(c['known_id'][nearest.argmin()]),closest_trail_px=round(nearest.min(),3)))
    save(ROOT/'locked_oct6_holdout.csv',audits)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['long','short','ee']:p.add_argument(k,type=Path)
    run(p.parse_args())
