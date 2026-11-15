"""Plot the frozen 3-GHz cutouts, optical images, and DESI [O II] line regions."""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from reproduce import HERE, INPUT, catalog, image_cutout

m = json.loads((HERE/'measurements.json').read_text())
spectra = pd.read_csv(INPUT/'desi.csv')
fig, axs = plt.subplots(2,3,figsize=(12.5,6.9),gridspec_kw={'width_ratios':[1.15,.93,1.26]})
for row,(radio,optical,spec) in zip(m,axs):
    tag=row['rgz'].split('_')[1][:7]
    obs=spectra[spectra.targetid==row['targetid']].iloc[0]
    image,head=image_cutout(INPUT/(tag+'_vlass.fits'))
    f=catalog(INPUT/(tag+'first.dat'))
    cosdec=np.cos(np.radians(obs.mean_fiber_dec))
    cx=(obs.mean_fiber_ra-float(head['CRVAL1']))*cosdec/float(head['CD1_1'])+float(head['CRPIX1'])-1
    cy=(obs.mean_fiber_dec-float(head['CRVAL2']))/float(head['CD2_2'])+float(head['CRPIX2'])-1
    rms=1.4826*np.median(np.abs(image-np.median(image)))
    radio.imshow(np.arcsinh(image/(2*rms)),origin='lower',cmap='gray_r',vmin=-.8,vmax=3.8,
                 extent=[cx+.5,cx+.5-image.shape[1],-.5-cy,image.shape[0]-.5-cy])
    xs=(f._RAJ2000-obs.mean_fiber_ra)*cosdec*3600
    ys=(f._DEJ2000-obs.mean_fiber_dec)*3600
    radio.scatter(xs,ys,s=78,facecolors='none',edgecolors='#cf5744',linewidths=1)
    radio.scatter([0],[0],marker='+',s=44,c='#2173aa')
    radio.set(xlim=(100,-100),ylim=(-65,65),xlabel='East offset (arcsec)',ylabel='North offset (arcsec)')
    radio.set_title(row['rgz'].replace('RGZ_','RGZ ')+'\nVLASS image · FIRST components',fontsize=10)
    optical.imshow(plt.imread(INPUT/(tag+'_optical.jpg')),extent=[13.1,-13.1,-13.1,13.1])
    core=f[f.FIRST==row['core']].iloc[0]
    dx=(core._RAJ2000-obs.mean_fiber_ra)*cosdec*3600
    dy=(core._DEJ2000-obs.mean_fiber_dec)*3600
    optical.scatter([dx],[dy],facecolors='none',edgecolors='#ef885b',s=78,label='FIRST core')
    optical.scatter([0],[0],marker='+',s=55,color='#8be7ff',label='DESI fiber')
    optical.set(xlim=(9,-9),ylim=(-9,9),xlabel='East offset (arcsec)',ylabel='North offset (arcsec)')
    optical.set_title('Legacy Surveys DR10',fontsize=10)
    optical.legend(fontsize=7,loc='lower left',framealpha=.8)
    d=pd.read_csv(INPUT/(tag+'_line.csv'));mu=3727.42*(1+obs.z);d=d[(abs(d.wavelength-mu)<47)&(d['mask']==0)&(d.ivar>0)]
    spec.step(d.wavelength,d.flux,where='mid',color='#96a1a8',linewidth=.8,label='DESI coadd')
    spec.plot(d.wavelength,gaussian_filter1d(d.flux,1.4),color='#245f79',linewidth=1.6,label='Smoothed')
    for line in [3727.092,3729.875]:
        spec.axvline(line*(1+obs.z),c='#c4573d',linestyle=':',linewidth=1)
    spec.set_title(f'[O II] at z = {obs.z:.5f}',fontsize=10)
    spec.set_xlabel('Observed vacuum wavelength (Å)')
    spec.set_ylabel(r'Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
    spec.legend(fontsize=7,frameon=False)
    for ax in [radio,optical,spec]:ax.tick_params(labelsize=8)
fig.tight_layout()
fig.savefig(HERE/'fields.png',dpi=155,bbox_inches='tight')
