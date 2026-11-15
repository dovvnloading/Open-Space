"""Make matched-scale 3-GHz image strips in mJy per beam."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

ROOT=Path(__file__).resolve().parent
with (ROOT/'observations.csv').open() as stream:
    records=list(csv.DictReader(stream))
names=['V1708 Sco','PGIR22akgylf','V1703 Sco']
fig, axes=plt.subplots(3,4,figsize=(9,6.7),layout='constrained')
for i,name in enumerate(names):
    rows=[r for r in records if r['name']==name]
    for ax in axes[i]:
        ax.set_axis_off()
    for j,r in enumerate(rows):
        ax=axes[i,j]
        path=ROOT/'cutouts'/r['cutout']
        with fits.open(path) as hdus:
            h=hdus[0].header
            data=hdus[0].data.squeeze()*1000
            w=WCS(h).celestial
        x,y=w.world_to_pixel_values(float(r['ra']),float(r['dec']))
        ix,iy=int(np.round(x)),int(np.round(y))
        image=ax.imshow(data[iy-12:iy+13,ix-12:ix+13],origin='lower',
                        cmap='inferno',vmin=-0.3,vmax=[2.8,3.4,1.2][i])
        # Plot the four sides of an open optical-position reticle.
        xx,yy=12+x-ix,12+y-iy
        for offset,finish in [(-4,-2.4),(2.4,4)]:
            ax.plot([xx+offset,xx+finish],[yy,yy],c='cyan',lw=0.9)
            ax.plot([xx,xx],[yy+offset,yy+finish],c='cyan',lw=0.9)
        ax.set_title(h['DATE-OBS'][:10],fontsize=9)
        if j==0:
            ax.text(-0.10,0.5,'AT 2022sfe' if name=='PGIR22akgylf' else name,
                    rotation=90,ha='right',va='center',transform=ax.transAxes)
    cb=fig.colorbar(image,ax=list(axes[i]),shrink=0.74)
    cb.set_label('mJy/beam',fontsize=8)
fig.savefig(ROOT/'radio-epochs.png',dpi=150)
