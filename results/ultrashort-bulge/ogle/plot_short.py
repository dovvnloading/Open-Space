"""Show actual peak-window sampling and competing models of four OGLE catalog alerts."""
import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import VBMicrolensing
from vet_short import HERE, lens


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=pathlib.Path); cache=ap.parse_args().cache
    vb=VBMicrolensing.VBMicrolensing()
    vb.LoadESPLTable(str(pathlib.Path(VBMicrolensing.__file__).parent/'data/ESPL.tbl'))
    plt.rcParams['svg.fonttype']='none'
    fig,axes=plt.subplots(2,2,figsize=(9,5.7),layout='constrained')
    windows=[(-.4,.4),(-.8,1.1),(-.34,.16),(-1.6,1.3)]
    for ax,rec,(left,right) in zip(axes.flat,json.loads((HERE/'short-fits.json').read_text()),windows):
        event=rec['event'];yr,_,n=event.split('-')
        t,m,err,_,_=lens.read_phot(cache/yr/f'blg-{n}.dat')
        t-=rec['assumed_t0_hjd']-2450000
        near=(t>=left)&(t<=right)
        ax.errorbar(t[near],m[near],yerr=rec['error_scale']*np.hypot(err[near],.004),fmt='o',ms=3,
                    color='#252525',lw=.75,label='OGLE I')
        grid=np.linspace(left,right,500)
        for k,color,style in [('FSPL','#177869','-'),('binary','#ba511d','--')]:
            f=rec['fits'][k];dt,lu,lt,*x=f['params']
            if k=='FSPL':a=vb.ESPLLightCurve([lu,lt,dt,x[0]],grid.tolist())[0]
            else:a=vb.BinaryLightCurve([x[0],x[1],np.exp(lu),x[2],x[3],lt,dt],grid.tolist())[0]
            model=lens.ZP-2.5*np.log10(np.array(a)*f['fs']+f['fb'])
            ax.plot(grid,model,style,color=color,lw=1.3,label='Finite-source point lens' if k=='FSPL' else 'Binary lens trial')
        ax.set_xlim(left,right)
        # Anchor the view to the data, since unsampled caustics can be arbitrarily bright.
        pad=max(float(np.ptp(m[near]))*.22,.035)
        ax.set_ylim(m[near].max()+pad,m[near].min()-pad)
        ax.set_title('OGLE-'+event,loc='left',fontsize=10,fontweight='bold')
        ax.set_xlabel('Days from EWS peak');ax.set_ylabel('I (mag)')
        ax.spines[['top','right']].set_visible(False)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside upper center',ncols=3,frameon=False)
    fig.savefig(HERE/'short-alerts.svg',metadata={'Date':None},bbox_inches='tight')

if __name__=='__main__':main()
