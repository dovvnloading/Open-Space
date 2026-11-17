"""Single-visit, two-filter Neptune orbit screen and Hippocamp forced control."""
import argparse
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from geometry import generate
from search import Cuts, file_hash, load_npz, prepare, sample_paths, snr, validate_cube, validate_tracks


def photometry(n,w):
    count=(w>0).sum(1)
    total,weight=n.sum(1),w.sum(1)
    flux=np.divide(total,weight,out=np.zeros_like(total),where=weight>0)
    drop=snr(total[:,None]-n,np.maximum(weight[:,None]-w,0))
    worst=np.min(np.where(w>0,drop,np.inf),axis=1)
    return count,snr(total,weight),np.where(count>1,worst,0),flux


def band_measure(data,tracks):
    cube,maps,cut=data
    validate_tracks(tracks,cube)
    return sample_paths(maps,tracks['x_path'],tracks['y_path'],cut)


def run(paths,out,grid,control):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    data=[]; cuts=Cuts(psf_support=.8)
    for p in paths:
        c=load_npz(p)
        if len(np.unique(c['visit']))!=1 or c['unmasked_control_body']!=str(control):
            raise ValueError('expected one-visit matched control or blind data')
        data.append((c,prepare(c,validate_cube(c,min_visits=1),cuts),cuts))
    if len({str(d[0]['visit'][0]) for d in data})!=1 or len({str(d[0]['band']) for d in data})!=2:
        raise ValueError('expected same visit in two different passbands')
    metadata={'cube_sha256':[file_hash(p) for p in paths],'bands':[str(d[0]['band']) for d in data],
              'visit':str(data[0][0]['visit'][0]),'psf_support':cuts.psf_support}
    if control is not None:
        rows=[]
        for c,maps,cut in data:
            path=c['known_path'][list(c['known_id']).index(control)]
            # Fixed on-sky offset tracks. Local comparisons are correlated.
            shifts=np.array([(0.,0.)]+[(r*np.cos(t),r*np.sin(t)) for r in [8.,12.,16.] for t in np.deg2rad(np.arange(0,360,15))])
            x=path[None,:,:,0]+shifts[:,0,None,None];y=path[None,:,:,1]+shifts[:,1,None,None]
            n,w=sample_paths(maps,x,y,cut);cnt,s,drop,flux=photometry(n,w)
            match=(cnt[1:]==cnt[0]);rows.append(dict(band=str(c['band']),n=int(cnt[0]),snr=float(s[0]),dropout=float(drop[0]),
                fit_uJy=float(flux[0]*1e6),n_controls=int(match.sum()),n_as_high=int((s[1:][match]>=s[0]).sum()),
                individual_scores=snr(n[0],w[0]).tolist(),exposure_id=c['exposure_id'].tolist()))
        (out/f'control_{control}.json').write_text(json.dumps({**metadata,'forced_body':control,'fits':rows},indent=2)+'\n')
        print(rows[0]['snr'],rows[1]['snr']);return
    lo,hi,step,dphase=grid;epoch=None
    with TemporaryDirectory() as work, (out/'radii.csv').open('w',newline='') as summary,(out/'passing.csv').open('w',newline='') as hits:
        radius_rows=csv.writer(summary);candidates=csv.writer(hits)
        radius_rows.writerow(['a_km','trials','eligible_both','pass_both','inverse_pass_both','max_joint_rank','filter1_ge_7','filter2_ge_7'])
        candidates.writerow(['a_km','phase_deg','filter1_score','filter2_score','dropout1','dropout2','n1','n2','f1_uJy','f2_uJy','joint_score'])
        for a in np.arange(lo,hi+step*1e-9,step):
            bands=[]
            for p,d in zip(paths,data):
                fname=Path(work)/'tracks.npz';generate(p,fname,grid=(a,a,step,dphase),epoch=epoch)
                t=load_npz(fname);epoch=float(t['epoch_jd_tdb'])
                n,w=band_measure(d,t);bands.append((photometry(n,w),photometry(-n,w)))
            def both(sign):
                c1,s1,drop1,_=bands[0][sign];c2,s2,drop2,_=bands[1][sign]
                eligible=(c1>=3)&(c2>=3)
                rank=np.minimum.reduce((s1/5,s2/5,drop1/3,drop2/3,np.hypot(s1,s2)/8))
                return eligible,np.where(eligible,rank,0)
            eligible,rank=both(0);_,neg=both(1)
            (n1,s1,d1,f1),(n2,s2,d2,f2)=bands[0][0],bands[1][0]
            radius_rows.writerow([a,len(rank),eligible.sum(),(rank>=1).sum(),(neg>=1).sum(),rank.max(),
                                  ((s1>=7)&(n1>=3)&(d1>=3)).sum(),((s2>=7)&(n2>=3)&(d2>=3)).sum()])
            for i in np.flatnonzero(rank>=1):
                candidates.writerow([a,t['phase_deg'][i],s1[i],s2[i],d1[i],d2[i],n1[i],n2[i],f1[i]*1e6,f2[i]*1e6,np.hypot(s1[i],s2[i])])
            summary.flush();hits.flush()
    (out/'parameters.json').write_text(json.dumps({**metadata,'epoch_jd_tdb':epoch,'grid':grid,
        'required_each_filter':{'exposure_count':3,'score':5,'worst_leave_one_out':3},'quadrature_min':8},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('f140m');p.add_argument('f210m');p.add_argument('output')
    p.add_argument('--grid',type=float,nargs=4,default=[30000,130000,500,1])
    p.add_argument('--control',type=int)
    a=p.parse_args();run([a.f140m,a.f210m],a.output,a.grid,a.control)
