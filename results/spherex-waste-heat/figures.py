"""Plot measured Hephaistos spectra and the WISE recovery surface."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np

ROOT = Path(__file__).resolve().parent

def rows(path):
    with (ROOT / path).open() as stream:
        return list(csv.DictReader(stream))

def array(data, key):
    return np.asarray([float(r.get(key) or np.nan) for r in data])

def plot_spectra():
    catalogs = rows('wise/catalog_photometry.csv')
    spec = rows('spectra/channel_spectra.csv')
    fig, grid = plt.subplots(5,2,figsize=(9.4,10.8),sharex=True)
    for target, ax in zip(catalogs, grid.T.flat):
        label = target['label']
        data = [r for r in spec if r['label']==label]
        wave, flux, err = [array(data,k) for k in ('wavelength_um','aperture_mJy','sky_sigma_mJy')]
        good = np.isfinite(flux) & np.isfinite(err)
        ax.vlines(wave[good],flux[good]-err[good],flux[good]+err[good],lw=.48,color='#719cb2')
        ax.scatter(wave,flux,s=5,color='#295e79',alpha=.86)
        for i,(w,zp) in enumerate(zip((3.4,4.6,12,22),(309.54,171.787,31.674,8.363)),1):
            mag = target[f'allwise_w{i}mpro']
            if not mag:
                continue
            point = 1000*zp*10**(-.4*float(mag))
            error = target[f'allwise_w{i}sigmpro']
            if target['allwise_ph_qual'][i-1] == 'U':
                ax.plot(w,point,'v',color='#aa5623',ms=4)
            elif error:
                ax.errorbar(w,point,yerr=point*.4*np.log(10)*float(error),fmt='D',ms=4,elinewidth=.9,capsize=2,color='#aa5623')
        ax.axhline(0,lw=.45,color='#90989e')
        ax.set(xscale='log',xlim=(.68,27))
        ax.set_yscale('symlog',linthresh=.45)
        ax.yaxis.set_major_locator(FixedLocator([-10,-1,0,1,10,100]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_: f'{x:g}'))
        ax.set_xticks([1,2,3,5,12,22],labels=['1','2','3','5','12','22'])
        ax.grid(axis='x',lw=.5,alpha=.2)
        ax.margins(y=.12)
        ax.text(.035,.06,'Hephaistos '+label,transform=ax.transAxes,ha='left',va='bottom',fontweight='bold',fontsize=9)
        ax.tick_params(labelsize=8)
    fig.suptitle('Infrared light at the ten Hephaistos positions',fontsize=17,y=.995)
    handles=[Line2D([],[],marker='.',color='#295e79',ls='none',label='SPHEREx DR1, strict mask'),
             Line2D([],[],marker='D',color='#aa5623',ls='none',ms=5,label='AllWISE catalog')]
    fig.legend(handles=handles,loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,.969))
    fig.supylabel('Flux density (mJy, symmetric log scale)',x=.015,fontsize=10)
    fig.supxlabel('Observed wavelength (µm)',y=.025,fontsize=10)
    fig.text(.5,.003,'DR1 bars show local blank-aperture scatter where measurable. AllWISE bars show catalog errors. Both beams include neighboring light.',ha='center',fontsize=8)
    fig.tight_layout(rect=(.035,.043,1,.942),h_pad=.4,w_pad=1)
    fig.savefig(ROOT/'target_spectra.png',dpi=130)
    plt.close(fig)

def plot_recovery():
    data = rows('blind/injections.csv')
    temps = sorted({int(r['dust_K']) for r in data})
    levels = sorted({float(r['L_excess_over_L_star']) for r in data})
    pct = {(int(r['dust_K']),float(r['L_excess_over_L_star'])):100*float(r['recovery_fraction']) for r in data}
    a = np.array([[pct[t,l] for l in levels] for t in temps])
    fig,ax=plt.subplots(figsize=(8,3.8))
    ax.imshow(a,vmin=0,vmax=100,aspect='auto',cmap='Blues',interpolation='nearest')
    for i,j in np.ndindex(a.shape):
        ax.text(j,i,f'{a[i,j]:.0f}%',ha='center',va='center',color='white' if a[i,j]>65 else '#16324b',fontsize=9)
    ax.set_xticks(range(len(levels)),[f'{l*100:g}%' for l in levels])
    ax.set_yticks(range(len(temps)),temps)
    ax.set_xlabel('Injected luminosity as a fraction of stellar luminosity')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('W3 or W4 recovery on 334 measured control stars',loc='left',pad=10,fontweight='bold',fontsize=13)
    fig.tight_layout()
    fig.savefig(ROOT/'wise_recovery.png',dpi=130)
    plt.close(fig)

if __name__=='__main__':
    plot_spectra()
    plot_recovery()
