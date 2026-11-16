"""Conditional scene-model sensitivity audit, NOT a calibrated molecular retrieval.

Fits each image independently; spectra are not tied across exposures. Neighbor
amplitudes are free and nonnegative, background is a plane. Radiometric ambiguity
and bright-source noise remain. See generated report for interpretation.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy import log
from scipy.ndimage import map_coordinates, maximum_filter
from scipy.signal import fftconvolve
from scipy.optimize import lsq_linear
from scipy.integrate import trapezoid
log.setLevel('ERROR')


def solve(columns, image, variance, mask, positive):
    A=np.column_stack([c[mask] for c in columns]); err=np.sqrt(variance[mask])
    fit=lsq_linear(A/err[:,None],image[mask]/err,bounds=([0.]*positive+[-np.inf]*(len(columns)-positive),np.full(len(columns),np.inf)),tol=1e-11)
    pred=sum(c*b for c,b in zip(columns,fit.x))
    return fit.x,pred,float(np.sum(((image[mask]-pred[mask])/err)**2))


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);p.add_argument('--calibration',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    rdir=a.root/'results/spherex-dr1-moving-spectra'; manifest=json.loads((rdir/'swan-full-d5.json').read_text())
    rows=[r for r in manifest if r['start'].startswith('2025-12-25') or r['start'].startswith('2026-01-02 10:')]
    bc=fits.getdata(a.calibration/'ParameterFiles/20260604_SSDC_BC_Band5.fits',1);ch=fits.getdata(a.calibration/'ChannelMaps/spectral_channels_17ch_2040.fits',1)
    summary=np.genfromtxt(next((a.calibration/'FidSpecSummary/SubChannel17/Band5').glob('*.csv')),delimiter=',',names=True)
    spec=np.loadtxt(rdir/'psg/co2_shape_100K_response.txt');wave=spec[:,0];fnu=wave**2*spec[:,-1]; fnu/=np.max(fnu)
    curves={k:np.genfromtxt(a.calibration/f'MiddlePixel/SubChannel17/Band5/FiducialSpec_Band5_SubCh{k:02d}_17_Cent.csv',delimiter=',',names=True) for k in range(1,18)}
    def response(bcent,channel):
        result=np.zeros_like(bcent)
        for k in np.unique(channel):
            if k<1 or k>17:continue
            m=channel==k; s=curves[k];grid=np.linspace(bcent[m].min()-.001,bcent[m].max()+.001,80);vals=[]
            for b in grid:
                w=s['Wavelength_um']+b-summary['BandCenter_um'][k-1];R=s['NormTransmission'];vals.append(trapezoid(R*np.interp(w,wave,fnu,left=0,right=0),w)/trapezoid(R,w))
            result[m]=np.interp(bcent[m],grid,vals)
        return result
    out=[]
    for row in rows:
        path=a.cache/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')
        if not path.exists():out.append({'start':row['start'],'error':'missing input'});continue
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        if actual!=row['sha256']:raise ValueError('Input hash mismatch '+str(path))
        with fits.open(path) as h:
            im=h['IMAGE'].data.astype(float);var=h['VARIANCE'].data;fl=h['FLAGS'].data; x,y=row['x'],row['y']; yy,xx=np.indices(im.shape);rad=np.hypot(xx-x,yy-y)
            good=np.isfinite(im)&np.isfinite(var)&(var>0)&((fl&~(2097152|524288))==0)
            w=WCS(h['IMAGE'].header,key='A');w.sip=None; ax,ay=w.wcs_pix2world(x,y,0);t=h['EPSF'].data;k=np.argmin((t['XCENTER']-ax)**2+(t['YCENTER']-ay)**2);ps=t['EPSF'][k].astype(float);ps/=ps.sum()
            def point(sx,sy):return 25*map_coordinates(ps,[(yy-sy)*5+16,(xx-sx)*5+16],order=3,mode='constant')
            neighbors=[];origin='image-selected peaks outside 3 pixels, not verified stars'
            if row['start'].startswith('2025-12-25'):
                tab=ascii.read(rdir/'allwise-dec25.tbl',format='ipac');sw=WCS(h['IMAGE'].header)
                nx,ny=sw.world_to_pixel_values(tab['ra'],tab['dec']);neighbors=[(float(u),float(v)) for u,v in zip(nx,ny) if np.hypot(u-x,v-y)<13];origin='AllWISE positions; free contemporaneous amplitudes, no catalog flux subtraction'
            else:
                local=(im==maximum_filter(im,size=5))&good&(rad>3)&(rad<12)&(im>np.nanmedian(im)+5*np.sqrt(var))
                ids=np.argwhere(local);ids=sorted(ids,key=lambda q:im[tuple(q)],reverse=True)[:5];neighbors=[(float(q[1]),float(q[0])) for q in ids]
            # Oversampled intrinsic scene; response multiplication occurs before PSF convolution.
            gy,gx=np.indices((201,201));dx=(gx-100)/5;dy=(gy-100)/5;X=x+dx;Y=y+dy;AX,AY=w.wcs_pix2world(X,Y,0)
            centers=map_coordinates(bc,[AY,2039-AX],order=1);channels=ch[np.clip(np.rint(AY).astype(int),0,2039),np.clip(np.rint(AX).astype(int),0,2039)];q=response(centers,channels);q/=max(float(q[100,100]),1e-15)
            def conv_model(power,trail,angle,weighted):
                cm=np.zeros_like(dx)
                for offset in np.linspace(-trail/2,trail/2,9):
                    cm+=((dx-offset*np.cos(angle))**2+(dy-offset*np.sin(angle))**2+.1**2)**(-power/2)/9
                if weighted:cm*=q
                cv=fftconvolve(cm,ps,mode='same');model=map_coordinates(cv,[(yy-y)*5+100,(xx-x)*5+100],order=3,mode='constant');return model/model[rad<=5].sum()
            entry={'start':row['start'],'url':row['url'],'sha256':actual,'neighbors':neighbors,'neighbor_origin':origin,'fits':[]}
            stars=[point(*n) for n in neighbors];plane=[np.ones_like(im),xx-x,yy-y]
            for power in [.5,1.,1.5]:
              for trail,angle in [(0,0),(.5,0),(.5,np.pi/2)]:
               for weighted in [False,True]:
                coma=conv_model(power,trail,angle,weighted)
                for radius in [5.,8.,12.]:
                 for deblend in [False,True]:
                  cols=[coma]+(stars if deblend else [])+plane;mask=good&(rad<=radius);beta,pred,chi=solve(cols,im,var,mask,len(cols)-3)
                  residuals={}
                  for ap in [2.5,3.5,5.]:
                    m=good&(rad<=ap);residuals[str(ap)]={'residual_sum':float(np.sum((im-pred)[m])),'formal_sigma':float(np.sqrt(var[m].sum())),'coma_sum':float((coma*beta[0])[m].sum()),'n':int(m.sum())}
                  entry['fits'].append({'power':power,'trail_pixels':trail,'trail_angle':angle,'response_before_convolution':weighted,'radius':radius,'deblend':deblend,'parameters':beta.tolist(),'chi2':chi,'n':int(mask.sum()),'k':len(cols),'apertures':residuals})
            out.append(entry)
            print(row['start'],len(neighbors),'neighbors',flush=True)
    a.output.write_text(json.dumps({'scope':'Conditional per-image nuisance-model sensitivity, not joint multi-epoch spectral retrieval','limitations':['All pixels and epochs already inspected; no independent holdout','OUTLIER pixels included; no upstream flag simulation','Eq22/NormTransmission=throughput interpretation conditional; 100K template not temperature measurement','Response-before-convolution uses one local finite ePSF; wavelength-dependent optics unknown','0.5 pixel trails along two image axes are sensitivity bounds, not ephemeris-derived trails','Jan neighbors selected on same image can absorb coma; Dec catalog positions may have motion/variability','Independent formal variances only; fit statistics are descriptive, not detection significance','Finite 20 pixel intrinsic coma; unit flux within radius5 normalization; not production rates','Per-image free amplitudes do not enforce a common molecular spectrum across images'], 'observations':out},indent=2)+'\n')

if __name__=='__main__':main()
