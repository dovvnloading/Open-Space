"""Render the ZTF source and the two neighboring comparison stars."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from ground_photometry import BASE, GAIA, STARS, read_star, wrap


def main():
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 5), sharex=True, layout="constrained")
    ref = Time(2455197.5+GAIA['reference_time'], format='jd', scale='tcb')
    p = (Time(ref.jd+1/GAIA['frequency'], format='jd',scale='tcb').tdb-ref.tdb).jd
    for axcol, name, title in zip(axes.T, STARS, ['KIC 8110754 · binary','KIC 8110757 · KOI target','KIC 8110758 · neighbor']):
        t, band, oid, mag=read_star(name)
        phase=((t-ref.tdb.jd)/p-GAIA['primary_phase']+.25)%1-.25
        outside=(abs(wrap(phase))>.165)&(abs(wrap(phase+GAIA['primary_phase']-GAIA['secondary_phase']))>.15)
        for o in np.unique(oid):
            m=oid==o
            mag[m]-=np.median(mag[m&outside])
        for ax, f, color in zip(axcol,['zg','zr'],['#318670','#ad4a44']):
            keep=band==f
            ax.scatter(phase[keep],mag[keep],s=4,alpha=.15,color=color,edgecolors='none',rasterized=True)
            bins=np.linspace(-.25,.75,51)
            vals=[]
            for lo,hi in zip(bins[:-1],bins[1:]):
                group=keep&(phase>=lo)&(phase<hi)
                vals.append(np.median(mag[group]) if group.sum()>=3 else np.nan)
            ax.plot((bins[1:]+bins[:-1])/2,vals,color=color,lw=2)
            ax.axhline(0,color='.6',lw=.6,zorder=0)
            ax.set_xlim(-.25,.75);ax.set_ylim(.86,-.20)
            ax.spines[['top','right']].set_visible(False)
        axcol[0].set_title(title,fontsize=11)
        axcol[1].set_xlabel('Orbital phase')
    for ax, name in zip(axes[:,0],['ZTF g','ZTF r']):
        ax.set_ylabel(name+'  ·  Δ magnitude')
    fig.savefig(BASE/'ztf_eclipses.png',dpi=160)
    plt.close(fig)


if __name__=='__main__':
    main()
