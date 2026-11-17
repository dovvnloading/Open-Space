"""Plot the conditional projected H-band M/L for the grade-B follow-up target."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'inference'))
from infer import chi,luminosity,mass_floor,read

def main():
    target='102158888_2659753895653871052'
    solar=np.genfromtxt(ROOT/'inference/sun_shifted_band_ab.csv',delimiter=',',names=True)
    m=next(r for r in read(ROOT/'inference/mass_light.csv') if r['id_str']==target and r['band']=='h')
    mag=float(m['mag_E']);theta=float(m['theta_eff_arcsec'])
    fig,ax=plt.subplots(figsize=(7.8,4.5),layout='constrained')
    for zs,color in [(2,'#126b7c'),(3,'#7c51a1'),(4,'#af5c2b')]:
        z=np.arange(.08,min(2.76,zs-.01),.007)
        ml=np.array([mass_floor(theta,zl)/(1-chi(zl)/chi(zs))/luminosity(mag,zl,solar,'h') for zl in z])
        ax.plot(z,ml,color=color,lw=2.1,label=f'Source z = {zs}')
    for model,z,xoffset,text in [('Euclid_PHZ',.17,.12,'Euclid PHZ'),('Legacy_DR9',.32,.17,'Legacy DR9')]:
        r=next(v for v in read(ROOT/'inference/conditional_ML_posterior.csv')
               if v['id_str']==target and v['band']=='h' and v['z_l_model']==model and float(v['zs'])==2)
        med,lo,hi=map(float,(r['ML_p50'],r['ML_p16'],r['ML_p84']))
        ax.errorbar(z,med,yerr=[[med-lo],[hi-med]],fmt='o',color='#172b35',capsize=4,zorder=4)
        ax.annotate(text,(z,med),xytext=(z+xoffset,med*1.06),fontsize=9)
    ax.axvline(2.4573,color='#666666',lw=1.2,ls='--')
    ax.text(2.425,36,'Euclid physical-\nparameter fit',fontsize=9,va='center',ha='right')
    ax.set(yscale='log',ylim=(1.5,240),xlim=(0,2.78),xlabel='Assumed deflector redshift',ylabel='Projected H-band mass-to-light ratio (solar units)')
    ax.yaxis.grid(True,which='major',color='#dddddd')
    ax.legend(loc='lower left',ncol=3,frameon=False,fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    fig.savefig(ROOT/'grade_b_geometry.png',dpi=155,facecolor='white')
    plt.close(fig)
if __name__=='__main__':main()
