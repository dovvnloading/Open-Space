"""Analyze fixed sky geometry in independent follow-up images and along-axis catalogue photometry."""
from pathlib import Path
import argparse,json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt
from archive_audit import convert
from wake_quality import extended_profile

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('images',type=Path);p.add_argument('results',type=Path)
    p.add_argument('out',type=Path);a=p.parse_args();images=a.images
    r=json.loads((a.results/'xdf_w_606_profile.json').read_text())['detections'][59]
    h606=fits.open(next(images.glob('*f606w_v1_sci.fits')))[0];whst=WCS(h606.header)
    xy=np.array([r['x0'],r['y0']]);u=np.array([r['x1'],r['y1']])-xy;length=np.linalg.norm(u);u/=length
    cat=fits.getdata(next(images.glob('*cat.fits')),1)
    pix=whst.all_world2pix(np.column_stack((cat['RA'],cat['DEC'])),0)
    along=(pix-xy)@u*.06;cross=(pix-xy)@[-u[1],u[0]]*.06
    members=(along>=-2)&(along<=length*.06+2)&(abs(cross)<1.25)
    coords=[]
    for i in np.where(members)[0]:
        coords.append(dict(id=int(cat['ID'][i]),along=float(along[i]),cross=float(cross[i]),**{
            k:float(cat[k][i]) for k in ['RA','DEC','Z_BPZ','ZMIN_BPZ','ZMAX_BPZ','ODDS_BPZ','SPECZ','GRISMZ']}))
    bands=[]
    for band in ['f606w','f775w','f814w','f850lp']:
        h=fits.open(next(images.glob(f'*{band}_v1_sci.fits')))[0];bands.append((band,h.data,WCS(h.header).celestial,.06))
    jpath=images/'jades_west_f150w_sci.fits'
    if jpath.exists():
        h=fits.open(jpath)[0];bands.append(('JADES F150W',h.data,WCS(h.header).celestial,.02999474))
    for f in sorted(images.glob('legacy_west*')):
        h=fits.open(f)[0];bands.append((h.header['VERSION']+' r',h.data[1],WCS(h.header).celestial,.262))
    fig,axes=plt.subplots(2,4,figsize=(10,5.2),squeeze=False);diag=[]
    for (b,im,w,scale),ax in zip(bands,axes.ravel()):
        row=convert(r,whst,w);cx,cy=round(row['x']),round(row['y'])
        rad=round(4.8/scale);x0,y0=cx-rad,cy-rad
        z=im[y0:cy+rad+1,x0:cx+rad+1];valid=np.isfinite(z)
        scale_obs=bool(valid.sum()>0 and np.nanstd(z)>0)
        local=dict(row)
        for k in ('x0','x1'):local[k]-=x0
        for k in ('y0','y1'):local[k]-=y0
        item=dict(band=b,nontrivial_coverage=scale_obs)
        if scale_obs:
            if scale<.1:item.update(extended_profile(z,valid,local))
            img=gaussian_filter(np.nan_to_num(z),max(.045/scale,.45));m=np.median(img);s=1.4826*np.median(abs(img-m))
            ax.imshow(np.arcsinh((img-m)/max(s,1e-14)),origin='lower',cmap='gray',vmin=-1.5,vmax=4.5,
                extent=(x0-.5,cx+rad+.5,y0-.5,cy+rad+.5))
            ax.scatter([row['x0'],row['x1']],[row['y0'],row['y1']],s=23,facecolors='none',edgecolors='gold',lw=.7)
        ax.set_title(b);ax.set_axis_off();diag.append(item)
    for ax in axes.ravel()[len(bands):]:ax.set_axis_off()
    fig.tight_layout();fig.savefig(a.out/'western_followup.jpg',dpi=110);plt.close(fig)
    (a.out/'western_followup.json').write_text(json.dumps(dict(geometry='XDF west F606W #59, fixed axis',along_from='endpoint 0',
        catalog_half_width_arcsec=1.25,beyond_end_arcsec=2,along_catalog=coords,independent_bands=diag),indent=2)+'\n')
if __name__=='__main__':main()
