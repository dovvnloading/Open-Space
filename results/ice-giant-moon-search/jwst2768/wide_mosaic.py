"""Prepare one target-registered, robustly standardized wide NIRCam map per physical visit."""
import argparse,csv,hashlib,json
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage,optimize
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import ARCSEC,AU,sky_xy
from prepare import inverse_sky
from ephemeris import interpolate


def center_of_limb(sci,cov,wcs,ra,dec,radius):
    x,y=wcs.all_world2pix(ra,dec,0);n=int(np.ceil(radius+16))
    x0,y0=int(np.rint(x))-n,int(np.rint(y))-n
    im=sci[y0:y0+2*n+1,x0:x0+2*n+1].astype(float);weight=cov[y0:y0+2*n+1,x0:x0+2*n+1]
    good=(weight>0)&np.isfinite(im)
    sm=ndimage.gaussian_filter(np.where(good,im,0),1)/np.maximum(ndimage.gaussian_filter(good.astype(float),1),1e-4)
    th=np.linspace(0,2*np.pi,360,endpoint=False);rad=np.arange(radius-7,radius+7,.15)
    px=x-x0+np.cos(th[:,None])*rad;py=y-y0+np.sin(th[:,None])*rad
    v=ndimage.map_coordinates(sm,[py,px],order=1);grad=np.gradient(v,rad,axis=1)
    j=grad.argmin(1);angle=np.arange(len(j));use=(j>2)&(j<len(rad)-3)&(grad[angle,j]<0)
    if use.sum()<180:raise ValueError('insufficient visible limb')
    rx=(px[angle,j]-(x-x0))[use];ry=(py[angle,j]-(y-y0))[use]
    def residual(par):
        dx,dy,r,c,s=par;t=np.arctan2(ry-dy,rx-dx)
        return np.hypot(rx-dx,ry-dy)-r-c*np.cos(2*t)-s*np.sin(2*t)
    fit=optimize.least_squares(residual,[0,0,radius,0,0],loss='soft_l1',f_scale=.5)
    if np.hypot(*fit.x[:2])>6 or abs(fit.x[2]-radius)>6:raise ValueError('implausible limb solution')
    return float(x+fit.x[0]),float(y+fit.x[1]),int(use.sum()),float(np.median(abs(fit.fun)))


def block_stats(raw,weight,b):
    ny,nx=[n//b for n in raw.shape]; sl=np.s_[:ny*b,:nx*b]
    v=raw[sl];ok=np.isfinite(v)&np.isfinite(weight[sl])&(weight[sl]>0)
    dims=(ny,b,nx,b);n=ok.reshape(dims).sum(axis=(1,3))
    img=np.where(ok,v,0).reshape(dims).sum(axis=(1,3))/np.maximum(n,1)
    return img.astype('f4'),n>=b*b-1


def robust_map(img,valid):
    # Field stars (and bright co-moving moons) occupy a small fraction of
    # each 48x48-pixel tile. Remove the locally smooth signal before
    # measuring noise from the remaining complete pixels.
    m=valid.astype('f4');cov=ndimage.gaussian_filter(m,9);bg=ndimage.gaussian_filter(np.where(valid,img,0),9)/np.maximum(cov,1e-5)
    ok=valid&(cov>.5);z=ndimage.gaussian_filter(np.where(ok,img-bg,0),.58)
    cov=ndimage.gaussian_filter(ok.astype('f4'),.58);ok&=cov>.97
    tile=48;ny,nx=z.shape;gh,gw=(ny+tile-1)//tile,(nx+tile-1)//tile
    med=np.full((gh,gw),np.nan);noise=med.copy()
    for j in range(gh):
        for i in range(gw):
            box=np.s_[j*tile:(j+1)*tile,i*tile:(i+1)*tile];s=z[box][ok[box]]
            if len(s)>=120:
                med[j,i]=np.median(s);noise[j,i]=1.4826*np.median(abs(s-med[j,i]))
    keep=np.isfinite(noise)&(noise>0)
    nearest=ndimage.distance_transform_edt(~keep,return_distances=False,return_indices=True)
    yy,xx=np.indices(z.shape,dtype='f4');xy=np.array([(yy+.5)/tile-.5,(xx+.5)/tile-.5])
    mu=ndimage.map_coordinates(med[tuple(nearest)],xy,order=1,mode='nearest')
    sigma=ndimage.map_coordinates(noise[tuple(nearest)],xy,order=1,mode='nearest')
    return np.where(ok,(z-mu)/sigma,0).astype('f4'),ok


def run(files,ephem,receipts,out,source):
    rows=json.loads(Path(receipts).read_text())['tile_inputs']; frames={}
    for q in rows:frames[q['exposure_id'].rsplit('_',2)[0]]=q
    out=Path(out);out.mkdir(parents=True,exist_ok=True);results=[]
    for file in sorted(map(Path,files)):
        with fits.open(file,memmap=True) as hd:
            head=hd[0].header;s=hd['SCI'];h=s.header;w=WCS(h, naxis=2);cov=hd['WHT'].data
            p0=head['EXPSTART'];pe=head['EXPEND']
            subs=sorted([r for r in frames.values() if p0-1e-6 <= r['mjd_start'] and r['mjd_end']<=pe+1e-6],key=lambda q:q['mjd_start'])
            if len(subs)!=5:raise ValueError(f'{file}: expected five calibrated exposure times, got {len(subs)}')
            mjd=np.array([[(r['mjd_start']+r['mjd_end'])/2] for r in subs]).ravel()
            _,_,ids,obs=interpolate(ephem,mjd)
            b=3; pix=np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600
            x,y,nlimb,rms=center_of_limb(s.data,cov,w,h['MT_AVRA'],h['MT_AVDEC'],obs[0,:,2].mean()/pix/2)
            ra,dec=w.all_pix2world(x,y,0)
            # Linearized east/north tangent vectors have subpixel error checked
            # over the full 50-arcsecond search field for each visit.
            M=[]
            for ee,nn in [(1,0),(0,1)]:
                wr,wd=inverse_sky(ra,dec,ee/ARCSEC,nn/ARCSEC);dx,dy=w.all_world2pix(wr,wd,0)
                M.append([(dx-x)/b,(dy-y)/b])
            M=np.array(M);cen=(np.array([x,y])-(b-1)/2)/b
            werr=0.
            for ee,nn in [(-50,-50),(-50,50),(50,50),(50,-50)]:
                wr,wd=inverse_sky(ra,dec,ee/ARCSEC,nn/ARCSEC);a1,a2=w.all_world2pix(wr,wd,0)
                p=np.array([a1,a2]);linear=(np.array([ee,nn])@M)*b+[x,y]
                werr=max(werr,np.linalg.norm(linear-p)/b)
            if werr>.05:raise ValueError('inadequate TAN approximation at field edge')
            img,valid=block_stats(s.data,cov,b)
        z,good=robust_map(img,valid)
        px=pix*b; ygrid,xgrid=np.ogrid[:z.shape[0],:z.shape[1]]
        # Bright, extended planets and moons have diffraction and saturation
        # wings beyond their angular diameters. These radii define exclusions.
        fixed={701:5.,702:5.,703:6.,704:6.,705:3.,715:1.5}
        cutmask=good.copy()
        for k,body in enumerate(ids):
            off=sky_xy(obs[k,:,0],obs[k,:,1],obs[0,:,0],obs[0,:,1],1.,(0.,0.))
            xy=np.stack((-off[:,0],off[:,1]),axis=1)@M+cen
            radius=12. if k==0 else max(fixed.get(int(body),.45),float(obs[k,:,2].max())/2+.4)
            for xx,yy in xy:
                if -radius/px<=xx<=z.shape[1]+radius/px and -radius/px<=yy<=z.shape[0]+radius/px:
                    mask=(xgrid-xx)**2+(ygrid-yy)**2 <= (radius/px)**2
                    cutmask[mask]=False
        name=subs[0]['exposure_id'].split('_')[0]
        np.savez_compressed(out/(name+'.npz'),z=z,good=cutmask,raw_good=good,mjd=mjd,path_ra=obs[0,:,0],
                            path_dec=obs[0,:,1],path_delta=obs[0,:,4],M=M,center=cen,pixel_scale=px,
                            source_filename=file.name,visit=name)
        results.append(dict(visit=name,file=file.name,science_shape=list(img.shape),bin_native=b,scale_arcsec=px,
                      planet_native_x=x,planet_native_y=y,n_limb_samples=nlimb,median_limb_residual_native_px=rms,
                      tan_linearization_max_output_px=float(werr),start_mjd=p0,end_mjd=pe,
                      contributing_detector_exposures=[q['exposure_id'].rsplit('_',2)[0] for q in subs]))
    Path(source).write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for q in ['eph','cal_receipt','out','receipt']:p.add_argument(q)
    p.add_argument('files',nargs='+');a=p.parse_args();run(a.files,a.eph,a.cal_receipt,a.out,a.receipt)
