"""Show the two independently registered IRAC2 visits on one celestial grid."""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.wcs import WCS
from measure import HERE, SOURCE, sample

rows = [x for x in json.loads((HERE/'measurements.json').read_text()) if x['band'] == 'I2']
east, north = np.meshgrid(np.linspace(-9,20,291), np.linspace(-8,13,211))
ra = SOURCE[0] + east / (3600*np.cos(np.deg2rad(SOURCE[1])))
dec = SOURCE[1] + north / 3600
fig, panels = plt.subplots(1, 2, figsize=(8.7, 4.3), sharex=True, sharey=True)
for ax, row in zip(panels, rows):
    img, hdr = fits.getdata(HERE/'cache'/row['filename'], header=True)
    w = WCS(hdr).celestial
    dx, dy = row['wcs_shift_pix']
    x, y = w.world_to_pixel_values(ra, dec)
    data = sample(img, (x+dx, y+dy))
    sky = np.nanmedian(data)
    scatter = 1.4826*np.nanmedian(abs(data-sky))
    stretched = np.arcsinh((data-sky)/(2*scatter))
    ax.imshow(stretched, extent=[-9.05,20.05,-8.05,13.05], origin='lower',
              cmap='gray_r', vmin=-.7, vmax=3.5)
    for pos, marker, color in [(row['fixed_sky']['xy_pix'], 'o', '#0078ac'),
                               (row['comoving']['xy_pix'], 'D', '#b74813'),
                               (row['host_xy_pix'], '+', '#792b84')]:
        a, b = w.pixel_to_world_values(pos[0]-dx, pos[1]-dy)
        e = (a-SOURCE[0])*3600*np.cos(np.deg2rad(SOURCE[1]))
        n = (b-SOURCE[1])*3600
        ax.plot(e, n, marker=marker, markersize=14, markeredgewidth=1.6,
                linestyle='', color=color, markerfacecolor='none')
    ax.set_title(row['date'][:10] + '  ·  IRAC 4.5 µm')
    ax.set_xlabel('East of JWST source (arcsec)')
    ax.set_xlim(20,-9)
panels[0].set_ylabel('North of JWST source (arcsec)')
legend = [Line2D([],[],marker=marker,markerfacecolor='none',linestyle='',markersize=9,
                 color=color,label=label,markeredgewidth=1.7) for marker,color,label in
          [('o','#0078ac','Fixed sky'),('D','#b74813','Co-moving prediction'),('+','#792b84','White dwarf')]]
fig.legend(handles=legend, loc='lower center', ncol=3, frameon=False)
fig.tight_layout(rect=(0,.105,1,1))
fig.savefig(HERE/'irac-motion.png', dpi=170, facecolor='white')
