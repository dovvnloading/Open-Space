"""Check every retained inner track against the local circular ring-light profile."""
import argparse,csv
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import generate
from search import Cuts,load_npz,prepare,sample_paths,score,validate_cube


def run(cube,initial,output):
    init=list(csv.DictReader(open(initial)))
    c=load_npz(cube);cut=Cuts(psf_support=.8);maps=prepare(c,validate_cube(c),cut)
    offset=np.array([0,-15,-12,-10,-8,-6,6,8,10,12,15])
    groups={}
    for row in init:groups.setdefault(float(row['a_km']),[]).append(float(row['phase_deg']))
    rows=[]
    with TemporaryDirectory() as tmp:
        el=Path(tmp)/'elements.csv';bank=Path(tmp)/'t.npz'
        for a,phases in groups.items():
            ph=np.mod(np.array(phases)[:,None]+offset,360)
            unique,inv=np.unique(ph,return_inverse=True)
            inv=inv.reshape(ph.shape)
            trials=np.zeros((len(unique),6));trials[:,0]=a;trials[:,-1]=unique
            np.savetxt(el,trials,delimiter=',',comments='',header='a_km,e,i_deg,node_deg,peri_deg,mean_deg')
            generate(cube,bank,elements=el,model='circular',epoch=2460347.42)
            t=load_npz(bank);n,w=sample_paths(maps,t['x_path'],t['y_path'],cut)
            flux=np.divide(n,w,out=np.full_like(n,np.nan),where=w>0)
            flank=flux[inv[:,1:]]; ctrl=np.nanmedian(flank,axis=1);ok=(np.isfinite(flank).sum(1)>=2)
            i=inv[:,0];nn=np.where(ok,n[i]-ctrl*w[i],0);ww=np.where(ok,w[i],0)
            fit=score(nn,ww,c['visit'],cut)
            for k,p in enumerate(phases):
                rows.append([a,p,float(fit['rank'][k]),fit['snr'][k],*fit['visit_snr'][k],*fit['n'][k]])
    with open(output,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['a_km','phase_deg','local_rank','pooled_control_score','visit1_score','visit2_score','visit3_score','n1','n2','n3']);w.writerows(rows)
    print(len(rows),sum(r[2]>=1 for r in rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for q in ['cube','initial','output']:p.add_argument(q)
    a=p.parse_args();run(a.cube,a.initial,a.output)
