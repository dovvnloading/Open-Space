"""Plot real-image recovery and convert its sampled 50% brackets to illustrative sizes."""
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from reduce import ROOT,DOMAINS,read_csv,write_csv
sys.path.insert(0,str(ROOT.parent))
from photometry import diameter

EE=.906   # WFC3/UVIS2 F845M, 10 native pixels. Mack et al. 2026, Table 4.
def main():
    sun=json.loads((ROOT.parent/'manifests/solar_calibration.json').read_text())['solar_AB_at_1au']['F845M']
    fig,axes=plt.subplots(1,2,figsize=(9.5,3.6),sharey=True,layout='constrained')
    rows=[]
    for ax,planet in zip(axes,['Uranus','Neptune']):
        ap=read_csv(planet.lower()+'/completeness.csv')
        settings=json.loads((ROOT/planet.lower()/'summary.json').read_text())
        for k,((lo,hi),color) in enumerate(zip(DOMAINS[planet],['#245987','#c15828'])):
            sub=[r for r in ap if int(r['zone'])==k]
            x=np.array([float(r['AB_mag']) for r in sub])+2.5*np.log10(EE)
            y=np.array([float(r['fraction']) for r in sub])
            low=np.array([float(r['lower_95']) for r in sub]); high=np.array([float(r['upper_95']) for r in sub])
            ax.fill_between(x,low,high,color=color,alpha=.14)
            ax.plot(x,y,'o-',color=color,label=f'{lo:,}–{hi:,} km',ms=4)
            before=np.flatnonzero(y>=.5)[-1]; after=before+1
            # A bracket is more faithful to the seven tested fluxes than a precise fitted crossing.
            size=lambda mag,p:float(diameter(mag,sun,settings['median_delta_au'],settings['median_sun_au'],p))
            rows.append(dict(planet=planet,a_min_km=lo,a_max_km=hi,bright_AB_total=round(x[before],3),
                         faint_AB_total=round(x[after],3),recovered_bright=sub[before]['dropout_recovered'],
                         recovered_faint=sub[after]['dropout_recovered'],trials_per_level=sub[before]['drawn'],
                         D_min_km_p012_at_faint=round(size(x[after],.12),1),D_max_km_p004_at_bright=round(size(x[before],.04),1)))
        ax.axhline(.5,ls=':',color='#555555',lw=1)
        ax.set(title=planet,xlabel='F845M total AB magnitude',ylim=(-.025,1.035))
        ax.legend(frameon=False,fontsize=8)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Fraction of all planted orbits recovered')
    fig.savefig(ROOT/'completeness.png',dpi=145)
    write_csv(ROOT/'depth_brackets.csv',rows)


if __name__=='__main__':main()
