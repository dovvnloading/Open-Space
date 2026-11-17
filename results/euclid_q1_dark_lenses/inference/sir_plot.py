"""Plot the independent, observed NISP line checks from IRSA combined spectra."""
import sys
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from infer import HERE

def run(cache):
    layouts=[('2655306065654633189',17040,17920,[(4862.68,'Hβ'),(4960.295,'[O III] 4960'),(5008.24,'[O III] 5008')],2.5387751,
              'Grade A   ·   emission at z = 2.539'),
             ('2659753895653871052',15700,16260,[(6564.61,'Hα')],1.43598294,
              'Grade B   ·   proposed Hα at z = 1.436')]
    fig,axs=plt.subplots(2,1,figsize=(8.5,5.4),layout='constrained')
    for ax,(key,lo,hi,lines,z,title) in zip(axs,layouts):
        tab=ET.parse(cache/(key+'.vot'))
        v=np.array([[float(c.text or 'nan')for c in row]for row in tab.findall('.//{*}TR')])
        w,f,mask,q,var,ndith,e=v.T
        ok=(mask.astype(int)%2==0)&(mask<64)&(ndith>=2)&(e>0)&(w>lo)&(w<hi)
        ax.fill_between(w[ok],(f-e)[ok],(f+e)[ok],color='#c6d9de',step='mid',label='1σ per bin')
        ax.step(w[ok],f[ok],where='mid',color='#164c5f',lw=1.45,label='Combined SIR spectrum')
        ax.set(xlim=(lo,hi));ax.set_title(title,loc='left',fontsize=11,fontweight='bold')
        for rest,label in lines:
            x=rest*(1+z);ax.axvline(x,color='#944719',lw=.85,ls='--')
            ax.annotate(label,(x,.965),xycoords=('data','axes fraction'),xytext=(3,0),textcoords='offset points',
                        fontsize=8,color='#6e340f',va='top')
        ax.tick_params(labelsize=8);ax.spines[['top','right']].set_visible(False)
    axs[1].set_xlabel('Observed wavelength (Å)',fontsize=10)
    fig.supylabel('Flux density  (10⁻¹⁶ erg s⁻¹ cm⁻² Å⁻¹)',fontsize=10)
    axs[0].legend(loc='upper left',fontsize=8,frameon=False,bbox_to_anchor=(0,-.11),ncol=2)
    fig.savefig(HERE/'nisp_spectra.png',dpi=145)
    plt.close(fig)

if __name__=='__main__':
    run(Path(sys.argv[1]))
