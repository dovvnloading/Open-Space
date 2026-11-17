"""Draw the IRAC2 residuals with the two trial source positions overlaid."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import SymLogNorm
import numpy as np
from astropy.io import fits
from measure_l788 import HERE,sample

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cache',type=Path);args=ap.parse_args()
    rows=json.loads((HERE/'forced-photometry.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(8.6,4.45),layout='constrained')
    handles=[Line2D([],[],marker=m,mec=c,mfc='none',color=c,linestyle='none',markersize=8,label=l)
             for m,c,l in [('o','#00dcbd','Fixed sky'),('s','#ef9c23','Fixed stellar offset'),('+','#cccccc','Host')]]
    for ax,aid in zip(axes,['57886720','65885696']):
        pair=[r for r in rows if r['method']=='reflected_linear' and r['aor']==aid];a,b=pair;cx,cy=a['host_x'],a['host_y']
        im=fits.getdata(next(args.cache.glob('*'+aid+'*maic.fits')))
        y,x=np.mgrid[round(cy)-25:round(cy)+26,round(cx)-25:round(cx)+26]
        residual=im[y,x]-sample(im,2*cx-x,2*cy-y)
        art=ax.imshow(residual,origin='lower',cmap='RdBu_r',norm=SymLogNorm(.08,vmin=-6,vmax=6),
                      extent=(x.min()-.5,x.max()+.5,y.min()-.5,y.max()+.5),interpolation='nearest')
        for r,s,color in [(a,'o','#00dcbd'),(b,'s','#ef9c23')]:
            ax.scatter(r['x'],r['y'],s=95,marker=s,edgecolor='black',facecolor='none',linewidth=2.5)
            ax.scatter(r['x'],r['y'],s=95,marker=s,edgecolor=color,facecolor='none',linewidth=1.4)
        ax.plot(cx,cy,'k+',ms=10);ax.plot(cx,cy,'+',color='#cccccc',ms=8)
        ax.set_title(f"{a['date'][:10]}  ·  AOR {aid}",fontsize=10);ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlabel('30.6 arcsec across · mosaic axes',fontsize=9)
    fig.colorbar(art,ax=axes,shrink=.76,label='180°-subtracted intensity (MJy/sr)')
    axes[0].legend(handles=handles,fontsize=7,loc='lower left',framealpha=.85)
    fig.suptitle('GJ 4274 at 4.5 µm with Spitzer/IRAC',fontsize=13)
    fig.savefig(HERE/'irac2-loci.png',dpi=155)
