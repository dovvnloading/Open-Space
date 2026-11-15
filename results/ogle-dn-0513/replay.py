"""Rebuild the OGLE cross-match and test for repeated brightenings. Python 3, numpy, scipy."""
import argparse
import hashlib
import io
import json
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

BASE = 'https://ftp.astrouw.edu.pl/ogle/ogle4/'
INPUTS = {
    'map_catalog': 'microlensing_maps/table3.dat',
    'dn_catalog': 'OCVS/CV/tab1.dat',
    'field_table': 'microlensing_maps/table7.dat',
    'map_photometry': 'microlensing_maps/phot.tar.gz',
    'dn_photometry': 'OCVS/CV/phot.tar.gz',
}
TARGET='BLG515.30.119116'
SECOND=2456040.84

def read_inputs(cache):
    cache.mkdir(parents=True,exist_ok=True)
    data={}
    for label,path in INPUTS.items():
        file=cache/label
        if not file.exists():
            with urllib.request.urlopen(BASE+path,timeout=120) as f:
                file.write_bytes(f.read())
        data[label]=file.read_bytes()
    return data

def coordinates(ra,dec):
    ra,dec=np.radians([ra,dec])
    return np.column_stack((np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)))

def sexagesimal(value):
    parts=np.abs(np.array(value.split(':'),dtype=float))
    return np.sign(-1 if value.startswith('-') else 1)*np.sum(parts/[1,60,3600])

def rows(data):
    return [l.split() for l in data.decode().splitlines() if l.strip() and not l.startswith('#')]

def curves(archive,wanted,offset):
    result={}
    with tarfile.open(fileobj=io.BytesIO(archive),mode='r|gz') as tf:
        for member in tf:
            name=Path(member.name).stem
            if name not in wanted or not member.isfile():continue
            arr=np.loadtxt(io.BytesIO(tf.extractfile(member).read()))
            arr[:,0]+=offset
            result[name]=arr
    return result

def triples(data,t0,tE):
    t,m,e=data.T
    nights,night=np.unique(np.floor(t),return_inverse=True)
    weights=e**-2
    w=np.bincount(night,weights=weights)
    mag=np.bincount(night,weights=m*weights)/w
    mid=np.bincount(night,weights=t*weights)/w
    err=w**-.5
    masksize=max(25,5*tE)
    out=[]
    for i in range(len(nights)-2):
        sl=slice(i,i+3)
        if max(np.diff(mid[sl]))>2.25 or abs(mid[i+1]-t0)<masksize:continue
        dist=abs(t-mid[i+1])
        ref=(dist>=12)&(dist<80)&(abs(t-t0)>masksize)
        if ref.sum()<30:continue
        signed=(np.median(m[ref])-mag[sl])/err[sl]
        out.append((min(signed),sum(signed)/np.sqrt(3),mid[i+1]))
    return sorted(out,reverse=True)

def measurements(data,center):
    t,m,e=data.T
    core=abs(t-center)<1.3
    ref=(abs(t-center)<80)&(abs(t-center)>12)
    base=float(np.median(m[ref]))
    excess=base-m[core]
    sigma=e[core]
    return dict(baseline_I=base,baseline_N=int(ref.sum()),epochs=data[core].tolist(),
                mean_brightening=float(np.average(excess,weights=sigma**-2)),
                measurement_error=float(sum(sigma**-2)**-.5),
                minimum_night_z=float(min(excess/sigma)),
                summed_z_over_sqrt_n=float(sum(excess/sigma)/np.sqrt(len(sigma))))

def main(cache):
    data=read_inputs(cache)
    mlrows=rows(data['map_catalog'])
    dnrows=[r for r in rows(data['dn_catalog']) if r[0].startswith('OGLE-BLG-DN')]
    mlxyz=coordinates(*np.array([[float(r[5]),float(r[6])] for r in mlrows]).T)
    dnxyz=coordinates(*np.array([[15*sexagesimal(r[1]),sexagesimal(r[2])] for r in dnrows]).T)
    tree=cKDTree(mlxyz)
    radius=2*np.sin(np.radians(2/3600)/2)
    near=tree.query_ball_point(dnxyz,radius)
    pairs=[(dnrows[i],mlrows[j],float(np.degrees(2*np.arcsin(np.linalg.norm(dnxyz[i]-mlxyz[j])/2))*3600)) for i,js in enumerate(near) for j in js]
    takes={r[0] for _,r,_ in pairs}|{r[0] for r in mlrows if r[1].startswith('BLG515.')}
    ml=curves(data['map_photometry'],takes,0)
    dn=curves(data['dn_photometry'],{r[0] for r,_,_ in pairs},2450000)
    def summarize_run(run):
        return dict(triples_tested=len(run),strongest_min_night_z=float(run[0][0]),summed_z_over_sqrt3=float(run[0][1]),middle_hjd=float(run[0][2])) if run else dict(triples_tested=0)
    overlaps=[]
    for dr,mr,sep in pairs:
        a,b=ml[mr[0]],dn[dr[0]]
        t0,te=map(float,mr[9:11])
        rec=dict(dn=dr[0],dn_field_id=dr[3]+'.'+dr[4],map_field_id=mr[0],separation_arcsec=sep,map_t0=t0,map_te=te,rate_weight=float(mr[-2]))
        for label,curve in [('dn',b),('map',a)]:
            t,m,e=curve.T
            peak=abs(t-t0)<max(20,te)
            brightest=curve[peak][np.argmin(m[peak])]
            rec[label+'_brightest_hjd']=float(brightest[0])
            rec[label+'_brightest_I']=float(brightest[1])
            rec[label+'_baseline_I']=float(np.median(m[abs(t-t0)>max(25,5*te)]))
            rec[label+'_search']=summarize_run(triples(curve,t0,te))
        overlaps.append(rec)
    dn_target=next(d for d,r,_ in pairs if r[0]==TARGET)[0]
    summary={}
    for label,arr in [('map',ml[TARGET]),('dn',dn[dn_target])]:
        summary[label]=measurements(arr,SECOND)
    controls=[]
    for mr in mlrows:
        if not mr[1].startswith('BLG515.'):continue
        arr=ml[mr[0]];t,m,e=arr.T
        t0,te=map(float,mr[9:11])
        peak=abs(t-SECOND)<1.3
        ref=(abs(t-SECOND)<80)&(abs(t-SECOND)>12)
        if sum(peak)!=3 or sum(ref)<30 or np.any(abs(t[peak]-t0)<max(25,5*te)):continue
        run=triples(arr,t0,te)
        if not run:continue
        z=(np.median(m[ref])-m[peak])/e[peak]
        controls.append(dict(source=mr[0],min_epoch_z=float(min(z)),summed_z_over_sqrt3=float(sum(z)/np.sqrt(3)),max_triplet_min_z=float(run[0][0]),triples=len(run)))
    comp=[r for r in controls if r['source']!=TARGET]
    dn_positions=np.array([[15*sexagesimal(r[1]),sexagesimal(r[2])] for r in dnrows])
    counts={}
    for shift in [-60,-40,-20,20,40,60]:
        shifted_ra=dn_positions[:,0]+shift/(3600*np.cos(np.radians(dn_positions[:,1])))
        shifted=coordinates(shifted_ra,dn_positions[:,1])
        counts[shift]=int(tree.query_ball_point(shifted,radius,return_length=True).sum())
    fieldrows=[r for r in mlrows if r[1].startswith('BLG515.')]
    published=next(r for r in rows(data['field_table']) if r[0]=='BLG515')
    weights=np.array([float(r[-2]) for r in fieldrows])
    times=np.array([float(r[10]) for r in fieldrows])
    candidate=next(r for r in fieldrows if r[0]==TARGET)
    w,candidate_te=float(candidate[-2]),float(candidate[10])
    field=dict(number=len(fieldrows),weight_sum=float(weights.sum()),
               reconstructed_mean_te=float(np.average(times,weights=weights)),
               published_mean_te=float(published[9]),published_rate=float(published[5]),
               published_rate_error=float(published[6]),published_tau=float(published[3]),
               published_tau_error=float(published[4]),
               single_exclusion_rate_fraction=w/sum(weights),
               single_exclusion_tau_fraction=w*candidate_te/sum(weights*times))
    info=dict(inputs={k:{'url':BASE+INPUTS[k],'sha256':hashlib.sha256(v).hexdigest()} for k,v in data.items()},map_rows=len(mlrows),field_sensitivity=field,bulge_dn_rows=len(dnrows),ra_shifts_arcsec=counts,pairs=overlaps,target=summary,field_controls=dict(number=len(comp),sum_triples=sum(r['triples'] for r in comp),max_min_z=max(r['max_triplet_min_z'] for r in comp),chip=[r for r in controls if r['source'].startswith('BLG515.30.')]),rate_weight_sum=sum(float(r[-2]) for r in mlrows))
    return info

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',default='ogle-downloads',type=Path)
    print(json.dumps(main(parser.parse_args().cache),indent=2,allow_nan=False))
