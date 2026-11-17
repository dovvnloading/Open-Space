"""Blind deep F150W2 trajectory search and independent F162M visit linkage."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import warnings

import numpy as np
root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root))
sys.path.insert(0,str(root/'jwst8975'))
import geometry
import search
from pilot_inputs import radial_psf
from ring_check import subtract_ring_profile


def remove_stationary(cube):
    data = np.where(cube['valid'],cube['image'],np.nan)
    residual = cube['image'].copy()
    valid=cube['valid'].copy()
    for i in range(len(data)):
        other = np.delete(data,i,axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            sky = np.nanmedian(other,axis=0)
        valid[i] &= np.isfinite(sky) & ((np.isfinite(other)).sum(axis=0)>=2)
        residual[i] -= np.nan_to_num(sky)
    return {**cube,'image':residual,'valid':valid}


def evaluate(cube, bank, psfs, maps, cuts):
    cache = search.TemplateCache(cube,bank,psfs)
    return cache.measure(maps,len(bank['x']),cuts)


def deep_score(num,info):
    half = len(num.T)//2
    sums, weights = num.sum(axis=1),info.sum(axis=1)
    z=search.snr(sums,weights)
    dropout = search.snr(sums[:,None]-num,np.maximum(weights[:,None]-info,0))
    drop = np.min(np.where(info>0,dropout,np.inf),axis=1)
    early = search.snr(num[:,:half].sum(axis=1),info[:,:half].sum(axis=1))
    late = search.snr(num[:,half:].sum(axis=1),info[:,half:].sum(axis=1))
    eligible=((info>0).sum(axis=1)>=3) & (info[:,:half].sum(axis=1)>0) & (info[:,half:].sum(axis=1)>0)
    rank=np.where(eligible,np.minimum.reduce([z/5,drop/3,early/2,late/2]),0)
    return rank,z,early,late,eligible


def make_cube(path, ee, deep):
    c=search.load_npz(path)
    band='F150W2' if deep else 'F162M'
    stamp=radial_psf(ee,float(c['scale_arcsec']),band)
    c.update(psf=stamp,psf_fraction=stamp.sum())
    return remove_stationary(c) if deep else subtract_ring_profile(c)


def trials(filename, cube_path, output, elements, epoch):
    np.savetxt(filename,elements,delimiter=',',comments='',header='a_km,e,i_deg,node_deg,peri_deg,mean_deg')
    geometry.generate(cube_path,output,elements=filename,model='circular',epoch=epoch)
    return search.load_npz(output)


def run(deep_path, confirm_path, ee_path, out, low, high, da, dphase):
    deep=make_cube(deep_path,ee_path,True)
    ref=make_cube(confirm_path,ee_path,False)
    kd=search.validate_cube(deep,min_visits=1)
    kr=search.validate_cube(ref)
    cuts=search.Cuts(psf_support=.8)
    md,mr=search.prepare(deep,kd,cuts),search.prepare(ref,kr,cuts)
    out.mkdir(parents=True,exist_ok=False)
    with TemporaryDirectory() as tmp, open(out/'seeds.csv','x',newline='') as fseeds, open(out/'refined.csv','x',newline='') as fref,open(out/'coverage.csv','x',newline='') as fcov:
        seed,cross,cov=csv.writer(fseeds),csv.writer(fref),csv.writer(fcov)
        seed.writerow(['a_km','phase_deg','deep_rank','deep_score','early_score','late_score'])
        cross.writerow(['a_km','phase_deg','deep_score','deep_early','deep_late','confirmation_visit','confirmation_score','confirmation_dropout','n_confirm','deep_flux_Jy','confirmation_flux_Jy'])
        cov.writerow(['rmin','rmax','n_coarse','deep_eligible','deep_positive','deep_negative','refinements','linked','ref_eligible'])
        temp=Path(tmp)
        radii=np.arange(low,high+da*1e-9,da)
        for start in range(0,len(radii),3):
            r=radii[start:start+3]
            geometry.generate(deep_path,temp/'deep.npz',grid=[r[0],r[-1],da,dphase])
            bank=search.load_npz(temp/'deep.npz')
            search.validate_tracks(bank,deep)
            epoch=float(bank['epoch_jd_tdb'])
            num,info=evaluate(deep,bank,kd,md,cuts)
            rank,z,early,late,elig=deep_score(num,info)
            neg=deep_score(-num,info)[0]
            peaks=np.flatnonzero(rank>=1)
            for k in peaks:
                seed.writerow([bank['a_km'][k],bank['phase_deg'][k],rank[k],z[k],early[k],late[k]])
            nlinks,eligible_ref,nref=0,0,0
            if len(peaks):
                # Cover time-baseline phase drift on the confirmation night.
                # Radius increments <= 25 km; retain native discovery epoch and geometry.
                dr=np.linspace(-da/2,da/2,int(np.ceil(da/25))+1)
                dp=np.array([-dphase/2,0,dphase/2])
                pairs=np.unique([(bank['a_km'][k]+offs,(bank['phase_deg'][k]+ph)%360)
                                 for k in peaks for offs in dr for ph in dp],axis=0)
                pairs=pairs[(pairs[:,0]>=low)&(pairs[:,0]<=high)]
                nref=len(pairs)
                for b in range(0,len(pairs),700):
                    block=pairs[b:b+700]; elements=np.zeros((len(block),6))
                    elements[:,0],elements[:,-1]=block.T
                    bd=trials(temp/'elem.csv',deep_path,temp/'rd.npz',elements,epoch)
                    br=trials(temp/'elem.csv',confirm_path,temp/'rr.npz',elements,epoch)
                    dn,di=evaluate(deep,bd,kd,md,cuts)
                    d_rank,d_snr,e,l,_=deep_score(dn,di)
                    rn,ri=evaluate(ref,br,kr,mr,cuts)
                    for visit in np.unique(ref['visit']):
                        use=ref['visit']==visit
                        ns,ws=rn[:,use],ri[:,use]
                        sm,sw=ns.sum(1),ws.sum(1)
                        n=(ws>0).sum(1)
                        eligible_ref+=int(((n>=2)&(d_rank>=1)).sum())
                        q=search.snr(sm,sw)
                        left=search.snr(sm[:,None]-ns,np.maximum(sw[:,None]-ws,0))
                        leave=np.min(np.where(ws>0,left,np.inf),axis=1)
                        keep=(d_rank>=1)&(n>=2)&(q>=5)&(leave>=2)
                        nlinks+=int(keep.sum())
                        for j in np.flatnonzero(keep):
                            cross.writerow([*block[j],d_snr[j],e[j],l[j],visit,q[j],leave[j],n[j],
                                            dn[j].sum()/di[j].sum(),sm[j]/sw[j]])
            cov.writerow([*r[[0,-1]],len(rank),elig.sum(),len(peaks),(neg>=1).sum(),nref,nlinks,eligible_ref])
            for fh in [fcov,fref,fseeds]:fh.flush()
    metadata=dict(deep_cube_sha256=search.file_hash(deep_path),confirmation_cube_sha256=search.file_hash(confirm_path),
                  ee_sha256=search.file_hash(ee_path),code_sha256=search.file_hash(Path(__file__)),epoch_tdb_jd=epoch,
                  rmin_km=low,rmax_km=high,coarse_radius_step_km=da,coarse_phase_step_deg=dphase,
                  max_confirmation_radius_step_km=25,psf_support=cuts.psf_support)
    (out/'parameters.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['deep','ref','ee','out']:
        p.add_argument(name,type=Path)
    p.add_argument('--grid',nargs=4,type=float,default=[54500,86000,250,.75])
    a=p.parse_args()
    run(a.deep,a.ref,a.ee,a.out,*a.grid)
