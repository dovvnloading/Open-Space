"""Fit two eclipse centers to TESS only. Compare with OGLE-only model predictions."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent
PHASE=(np.arange(120)+.5)/120
GROUPS={'2018_2019':[1,2,13],'2020':[27,28],'2023':[67,68],'2025':[94,95,96]}
FIELDS=['primary_phase','secondary_phase','separation','primary_depth_e_s','secondary_depth_e_s','primary_width','secondary_width']

def wrap(q):return (q+.5)%1-.5

def fit(S,N,wide=False):
    # Each sector supplies its native-cadence phase means; pool sectors by covered
    # orbital-cycle count, avoiding overweighting high-cadence sectors.
    A=[];W=[]
    for s,n in zip(S,N):
        ss=s.sum(axis=0);nn=n.sum(axis=0)
        arr=np.divide(ss,nn,out=np.zeros_like(ss),where=nn>0)
        weight=nn/max(nn.mean(),1)*len(s)
        A.append(arr);W.append(weight)
    w=np.sum(W,axis=0); y=np.average(np.array(A),axis=0,weights=np.array(W))
    def model(z):
        c1,c2,d1,d2,b,a,h=z[:7];w1,w2=z[7:] if wide else (.025,.025)
        return b+a*np.cos(2*np.pi*PHASE)+h*np.sin(2*np.pi*PHASE)-d1*np.exp(-.5*(wrap(PHASE-c1)/w1)**2)-d2*np.exp(-.5*(wrap(PHASE-c2)/w2)**2)
    start=[0,.53,1.2,1.2,0,0,0];lo=[-.09,.41,0,0,-10,-10,-10];hi=[.11,.61,25,25,10,10,10]
    if wide:start.extend([.025,.025]);lo.extend([.010,.010]);hi.extend([.060,.060])
    res=least_squares(lambda z:(model(z)-y)*np.sqrt(w/w.mean()),start,bounds=(lo,hi),max_nfev=120)
    c1,c2,d1,d2=res.x[:4];return np.r_[c1,c2,c2-c1,d1,d2,res.x[7:] if wide else [.025,.025]]

def pfield(d,m,t):
    mod=d['models'][m];th=mod['theta'];U=mod.get('best_period_yr');x=(t-2456000.)/365.25
    base=th[0]+th[1]*x/10.
    if m=='apsidal':disp=th[2]*np.cos(2*np.pi*x/U)+th[3]*np.sin(2*np.pi*x/U)
    else:
        n={'constant':1,'linear':2,'quadratic':3,'cubic':4}[m]
        disp=sum(th[2+j]*(x/10.)**j for j in range(n))
    return np.array([base-disp/2,.5+base+disp/2,.5+disp])

def process(draws):
    rng=np.random.default_rng(2106275)
    meta=json.loads((HERE/'tess_inputs.json').read_text())
    model=json.loads((HERE.parent/'orbit'/'orbit.json').read_text())
    out={'fields':FIELDS,'draws':draws,'seed':2106275,'groups':{}}
    with np.load(HERE.parent/'joint'/'inputs'/'orbital_folds.npz') as data:
        for name,sects in GROUPS.items():
            S=[];N=[]
            for s in sects:
                cnt=data[f's{s}_counts']
                # Use complete orbits covering >=65% of typical expected
                # cadences and >=55 phase bins, retaining eclipse opportunities.
                dt=meta[str(s)]['cadence_days']
                ok=(cnt.sum(axis=1)*dt>.65*1.1696719)&((cnt>0).sum(axis=1)>=min(55,int(.65*1.1696719/dt)))
                N.append(cnt[ok]);S.append(data[f's{s}_sums'][ok])
            byap={}
            for k,ap in enumerate(['3x3','four','one']):
                a=[s[:,:,k] for s in S]
                point=fit(a,N)
                byap[ap]={'point':point.tolist(),'free_width_point':fit(a,N,True).tolist()}
                if k==0:
                    sims=[]
                    for _ in range(draws):
                        idx=[rng.integers(0,len(t),len(t)) for t in a]
                        sims.append(fit([t[i] for t,i in zip(a,idx)],[n[i] for n,i in zip(N,idx)]))
                    byap[ap]['ci95']=np.percentile(sims,[2.5,97.5],axis=0).T.tolist()
            mid=np.average([np.mean(meta[str(s)]['range_btjd'])+2457000 for s in sects],weights=list(map(len,N)))
            preds={k:pfield(model,k,mid).tolist() for k in ['constant','linear','quadratic','cubic','apsidal']}
            out['groups'][name]={'n_orbits':sum(map(len,N)), 'orbit_count_by_sector':dict(zip(sects,map(len,N))),'hjd_mid_approx':mid,'apertures':byap,
                                'ogle_only_predicted_centers_and_separation':preds}
            print(name,sum(map(len,N)),np.round(byap['3x3']['point'][:3],4),np.round(byap['3x3']['ci95'][2],4),'alt',np.round([byap[k]['point'][2] for k in byap],4),'aps,cubic',np.round([preds[k][2] for k in ['apsidal','cubic']],4),flush=True)
    (HERE/'tess_timing.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':
    cli=argparse.ArgumentParser();cli.add_argument('--draws',type=int,default=600)
    process(cli.parse_args().draws)
