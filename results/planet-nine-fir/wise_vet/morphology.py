"""Inspect registered native W2 science and opposite-season subtraction shapes."""
import argparse,csv,json,sys,warnings
from pathlib import Path
import numpy as np
from scipy.ndimage import maximum_filter

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--wise',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parent);args=ap.parse_args()
    sys.path[:0]=[str(args.wise),str(args.wise.parent/'orbits'),str(args.wise.parent/'other')]
    from search import Frame,tangent
    rows=list(csv.DictReader((args.out/'individual.csv').open()))
    from rr_w2 import load_ledger
    rec=load_ledger(args.wise,'RR-h270')
    fs={x['id']:Frame(x,args.cache/'RR-h270') for x in rec['products'] if 'int_sha256' in x and 'msk_sha256' in x}
    axis=np.arange(-24,25,1.);xx,yy=np.meshgrid(axis,axis);rr=np.hypot(xx,yy)
    # North and east orientation; 1 arcsecond grid oversamples the native 2.75-arcsec pixels.
    raw=[];panels=[]
    for k in [635,727]:
        for visit,ref in [(0,1),(1,0),(2,1)]:
            regs=[];scis=[]
            for r in [r for r in rows if int(r['template'])==k and int(r['visit'])==visit]:
                track=(float(r['ra']),float(r['dec']))
                xy=tangent(track,xx.ravel(),yy.ravel())
                f=fs[r['frame']]
                sci=f.sample(xy)
                refs=np.array([g.sample(xy) for g in fs.values() if g.visit==ref])
                static=np.nanmedian(refs,axis=0);static[np.isfinite(refs).sum(axis=0)<2]=np.nan
                regs.append(sci-static);scis.append(sci)
            A=[]
            for seq in [scis,regs]:
                seq=np.array(seq);a=np.nanmedian(seq,axis=0);a[np.isfinite(seq).sum(axis=0)<2]=np.nan
                A.append(a.reshape(xx.shape))
            science,difference=A
            info=dict(template=k,visit=visit)
            for label,im in zip(['science','difference'],A):
                filled=np.where(np.isfinite(im),im,-np.inf)
                peaks=(filled==maximum_filter(filled,size=5))&(rr<20)&np.isfinite(im)&(im>0.)
                radii=rr[peaks]
                info[label+'_distance_to_nearest_positive_peak_as']=float(radii.min()) if len(radii) else None
                info[label+'_central_filtered']=float(im[24,24]) if np.isfinite(im[24,24]) else None
                core=(rr<=6.4)&np.isfinite(im);side=(rr>=8)&(rr<=11)&np.isfinite(im)
                info[label+'_peak_core']=float(im[core].max()) if np.any(core) else None
                info[label+'_ring_med']=float(np.median(im[side])) if np.any(side) else None
            panels.append((science,difference));raw.append(info)
    (args.out/'morphology.json').write_text(json.dumps(raw,indent=2)+'\n')
    from matplotlib import pyplot as plt
    from matplotlib.patches import Circle
    fig,axes=plt.subplots(4,3,figsize=(8,9),layout='constrained')
    for i,(science,diff) in enumerate(panels):
        row=(i//3)*2;col=i%3
        for j,img in enumerate([science,diff]):
            ax=axes[row+j,col];lim=np.nanpercentile(abs(img),96)
            ax.imshow(img,origin='lower',extent=(-24.5,24.5,-24.5,24.5),vmin=-lim,vmax=lim,cmap='RdBu_r')
            ax.add_patch(Circle((0,0),3.2,fc='none',ec='lime',lw=1))
            if col==0:ax.set_ylabel(f"ID {[635,727][i//3]} · "+['science','difference'][j])
            ax.set(xticks=[-20,0,20],yticks=[-20,0,20],aspect='equal')
            if row+j==0:ax.set_title(['June 2023','Dec 2023','June 2024'][col])
    fig.supxlabel('East offset (arcsec)');fig.supylabel('North offset (arcsec)')
    fig.savefig(args.out/'rr_w2_profiles.png',dpi=130)
    print(raw)

if __name__=='__main__':
    warnings.filterwarnings('ignore',message='All-NaN slice encountered')
    main()
