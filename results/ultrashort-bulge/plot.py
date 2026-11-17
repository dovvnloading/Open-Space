"""Display the known terrestrial-mass control and the separate KMT-0456 binary test."""
import argparse
import gzip
import json
import tarfile
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import fit

HERE=Path(__file__).resolve().parent


def lc(kind,v,ref,t,vb):
    t=np.asarray(t)-ref
    dt,lte,lu=v[:3]
    if kind=='fs':
        return np.asarray(vb.ESPLLightCurve([lu,lte,dt,v[3]],t.tolist())[0])-1
    return np.asarray(vb.BinaryLightCurve([v[4],v[5],np.exp(lu),v[6],v[3],lte,dt],t.tolist())[0])-1


def main(cache,out):
    vb=fit.engine()
    fig,axes=plt.subplots(1,3,figsize=(13.2,4.2),constrained_layout=True)
    ref=2457557.8
    result=json.loads((HERE/'controls/earth-recovery.json').read_text())
    vals=[]
    with gzip.open(HERE/'controls/earth-figure1.tsv.gz','rt') as f:
        for line in f:
            r=line.strip().split('\t')
            if len(r)>=5 and r[0].strip().isdigit(): vals.append((float(r[1]),float(r[2]),float(r[3]),r[4]))
    t,m,e,s=map(np.asarray,zip(*vals))
    ax=axes[0]
    for site,color in [('OGLE','#242a5a'),('KMTC04','#c85b15')]:
        nn=(s==site)&(abs(t-ref)<.2)
        flux,err=fit.lens.to_flux(m,e)
        base=np.median(flux[(s==site)&(abs(t-ref)>1)])
        ax.errorbar((t[nn]-ref)*24,flux[nn]/base-1,err[nn]/base,fmt='.',color=color,alpha=.75,label=site,ms=4)
    tt=np.linspace(-.19,.19,400)
    mdl=result['finite_no_blend']
    amp=vb.ESPLLightCurve([np.log(mdl['u0']),np.log(mdl['te_d']),mdl['t0']-ref,np.log(mdl['rho'])],tt.tolist())[0]
    ax.plot(tt*24,np.array(amp)-1,'-',color='#188268',label='Finite source',lw=2)
    ax.set(xlabel='Hours from HJD 2457557.8',ylabel='Excess flux / baseline flux',title='OGLE-2016-BLG-1928')
    ax.legend(fontsize=8,frameon=False)
    res=json.loads((HERE/'peaks/KMT-2024-BLG-0456.json').read_text());ref=res['reference_jd']
    with tarfile.open(cache/'staging/kmt/KMT-2024-BLG-0456.tar.gz') as tar:
        group=0
        for member in tar:
            if not member.name.endswith('_I.pysis'):continue
            dat=np.loadtxt(tar.extractfile(member),ndmin=2)
            if not len(dat):continue
            j=group;group+=1
            if member.name not in ('KMTS42_I.pysis','KMTC42_I.pysis'): continue
            site=member.name[:4]
            ax=axes[1 if site=='KMTS' else 2]
            time,df,ef=dat[:,:3].T
            valid=(ef>0)&np.isfinite(df)&np.isfinite(ef)
            off=valid&(abs(time-ref)>12)
            base=np.median(-df[off]); total=res['caps'][j]/res['cap_factor']
            near=valid&(abs(time-2460404.6)<1.15)&(ef<.5*total)
            ax.errorbar(time[near]-2460404,(-df[near]-base)/total,ef[near]/total,color='#222b44',fmt='.',ms=4,label='I-band data')
            grid=np.linspace(2460403.45,2460405.75,700)
            for k,col,label in [('fs','#188268','Single lens, finite source'),('bin','#bb5924','Binary lens')]:
                model=res['models'][k]
                ax.plot(grid-2460404,lc(k,model['vector'],ref,grid,vb)*model['fs'][j]/total,color=col,lw=1.6,label=label)
            ax.set(xlabel='HJD − 2460404 (days)',title='KMT-2024-BLG-0456 · '+('South Africa' if site=='KMTS' else 'Chile'),ylim=(-.6,3.9))
            if site=='KMTS':ax.legend(fontsize=8,frameon=False)
    for ax in axes:
        ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    fig.savefig(out,dpi=150);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--out',type=Path,default=HERE/'control-and-alternative.png');a=p.parse_args();main(a.cache,a.out)
