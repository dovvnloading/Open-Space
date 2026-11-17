"""Cross-match independent HST detections, measure longitudinal support, and render audit stamps."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates,gaussian_filter
from matplotlib import pyplot as plt
from wake_filter import compatible

def support(image, row):
    """Twelve axial bins, on-axis median minus flanking strips, locally scaled by parallel apertures."""
    x,y,L,w=row['x'],row['y'],row['length_px'],row['fwhm_px']
    t=np.deg2rad(row['pa_image_deg']);u=np.array([np.cos(t),np.sin(t)]);v=np.array([-u[1],u[0]])
    a=np.linspace(-.47*L,.47*L,12*8)
    def profile(shift):
        def strip(c):
            a2,c2=np.meshgrid(a,np.linspace(-.3*w,.3*w,5)+c+shift)
            xy=np.array([x,y])[:,None,None]+u[:,None,None]*a2+v[:,None,None]*c2
            return np.nanmedian(map_coordinates(image,xy[::-1],order=1).reshape(5,12,8),axis=(0,2))
        return strip(0)-.5*(strip(-1.8*w)+strip(1.8*w))
    on=profile(0)
    side=np.array([profile(k*w) for k in [-9,-7,-5,5,7,9]])
    noise=1.4826*np.nanmedian(np.abs(side-np.nanmedian(side)))
    return dict(n_positive=int((on>0).sum()),n_above_2scatter=int((on>2*noise).sum()),
                n_bins=len(on), contrast=on.tolist(),scatter=float(noise))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('images',type=Path);p.add_argument('results',type=Path)
    args=p.parse_args(); data=args.images;out=args.results
    bands=['606','814']
    sources=[fits.open(next(data.glob('*f'+band+'w_v1_sci.fits')),memmap=True) for band in bands]
    images=[h[0].data for h in sources];wcs=WCS(sources[0][0].header)
    a,b=[json.loads((out/('xdf_'+k+'.json')).read_text())['detections'] for k in bands]
    matches=[(i,j) for i,x in enumerate(a) for j,y in enumerate(b) if compatible(x,y) and compatible(y,x)]
    fig,axs=plt.subplots(int(np.ceil(len(matches)/3)),3,figsize=(11.5,3.3*int(np.ceil(len(matches)/3))),squeeze=False)
    results=[]
    for (i,j),ax in zip(matches,axs.ravel()):
        x,y=a[i],b[j];X,Y=int(x['x']),int(x['y'])
        radius=round(x['length_px']*.68)+10
        sl=np.s_[Y-radius:Y+radius+1,X-radius:X+radius+1]
        z=[im[sl] for im in images]
        z=[np.arcsinh(np.maximum(im-np.median(im),0)/(2*1.4826*np.median(abs(im-np.median(im))))) for im in z]
        rgb=np.dstack([z[1],(z[0]+z[1])/2,z[0]])
        ax.imshow(np.clip(rgb/max(2.5,np.percentile(rgb,99.9)),0,1),origin='lower',extent=(X-radius-.5,X+radius+.5,Y-radius-.5,Y+radius+.5))
        # Off-axis parallel line and open endpoint markers leave putative feature uncovered.
        t=np.radians(x['pa_image_deg']);sx=-np.sin(t)*6;sy=np.cos(t)*6
        ax.plot([x['x0']+sx,x['x1']+sx],[x['y0']+sy,x['y1']+sy],lw=.7,color='lime')
        ax.scatter([x['x0'],x['x1']],[x['y0'],x['y1']],s=30,facecolors='none',edgecolors='yellow',lw=.7)
        ax.set_title(f"F606W #{i} / F814W #{j}",fontsize=9); ax.set_axis_off()
        results.append(dict(f606_id=i,f814_id=j,radec=wcs.all_pix2world([[x['x'],x['y']]],0)[0].tolist(),
            endpoints_radec=wcs.all_pix2world([[x['x0'],x['y0']],[x['x1'],x['y1']]],0).tolist(),
            f606=support(images[0],x),f814=support(images[1],y)))
    for ax in axs.ravel()[len(matches):]:ax.set_axis_off()
    fig.tight_layout();fig.savefig(out/'xdf_matches.jpg',dpi=110);plt.close(fig)
    (out/'xdf_matches.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
