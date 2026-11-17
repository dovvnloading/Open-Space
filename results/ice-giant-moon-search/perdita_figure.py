"""Plot measured and URA184 Perdita tracks relative to Uranus in two visits."""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent

def main(args):
    rows=list(csv.DictReader((ROOT/'jwst6379_9482/psf_astrometry/astrometry.csv').open()))
    rows+=list(csv.DictReader((ROOT/'jwst9482/integration_check/displaced_astrometry.csv').open()))
    fig,axes=plt.subplots(1,2,figsize=(9,4.1),layout='constrained')
    for ax,path,label in zip(axes,[args.feb,args.oct],['February 2, 2025','October 7, 2025']):
        with np.load(path) as c:
            frames=c['exposure_id'].tolist();scale=float(c['scale_arcsec']);cx,cy=c['center_xy']
            body=np.flatnonzero(c['known_id']==725)[0];tracks=c['known_path'][body]
            texp=c['path_mjd']
        obs=[];pred=[]
        for r in rows:
            if r['frame'] not in frames:continue
            j=frames.index(r['frame']);t=float(r['mjd_utc'])
            obs.append([-(float(r['x_sky'])-cx)*scale,(float(r['y_sky'])-cy)*scale])
            px,py=[np.interp(t,texp[j],col) for col in tracks[j].T];pred.append([-(px-cx)*scale,(py-cy)*scale])
        obs,pred=map(np.array,[obs,pred])
        ax.plot(*pred.T,color='#b25425',lw=1.9,ls='--')
        ax.scatter(*obs.T,c='#164e71',s=9,linewidths=0)
        for k in np.unique(np.linspace(0,len(obs)-1,8,dtype=int)):
            ax.plot([pred[k,0],obs[k,0]],[pred[k,1],obs[k,1]],c='#9a9a9a',lw=.8,zorder=0)
        ax.set_aspect('equal',adjustable='datalim');ax.set_title(label)
        ax.set_xlabel('East of Uranus (arcsec)');ax.set_ylabel('North of Uranus (arcsec)')
        ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.13)
        med=np.median(np.linalg.norm(pred-obs,axis=1))
        ax.text(.03,.04,f'{len(obs)} positions   |   median separation {med:.2f}″',transform=ax.transAxes,fontsize=9,
                bbox={'facecolor':'white','edgecolor':'none','alpha':.8})
    fig.legend([Line2D([],[],marker='o',linestyle='',color='#164e71'),Line2D([],[],ls='--',color='#b25425')],
               ['Measured centroid','JPL URA184 prediction'],loc='outside upper center',ncol=2,frameon=False)
    fig.savefig(ROOT/'perdita_tracks.png',dpi=150)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['feb','oct']:p.add_argument(k,type=Path)
    main(p.parse_args())
