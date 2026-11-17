"""Native trailed-PSF photometry of the measured Perdita track in both deep visits."""
import argparse,csv,json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import warnings
from astropy.io import fits
import numpy as np
from scipy.optimize import least_squares
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'jwst9482'));sys.path.insert(0,str(base/'jwst8975'))
from link_visits import search,trials
from native_plants import Detector
import positive_control as pc


def fit_flux(data,kernel):
    yy,xx=np.indices(data.shape);yy=yy-(data.shape[0]-1)/2;xx=xx-(data.shape[1]-1)/2
    ok=np.isfinite(data)&np.isfinite(kernel)
    if ok.sum()<30 or (kernel*ok).sum()<.25:return None
    bg=ok&(np.hypot(xx,yy)>5)
    if bg.sum()<12:return None
    m=np.median(data[bg]);noise=max(1e-4,1.4826*np.median(abs(data[bg]-m)))
    def loss(par):return ((par[0]*kernel+par[1]+par[2]*xx+par[3]*yy-data)/noise)[ok]
    sol=least_squares(loss,[.4,m,0,0],loss='soft_l1',f_scale=2)
    return sol.x[0],ok.sum(),(kernel*ok).sum()


def visit(cube,map_path,centroids,meta,pixels,files,ee,flat):
    c=search.load_npz(cube);q=search.load_npz(map_path)
    records=json.loads(Path(str(cube)+'.inputs.json').read_text())
    rows=list(csv.DictReader(centroids.open()));good={(int(r['j']),int(r['k'])):r for r in rows if r['edge_parameter']=='0' and float(r['model_peak_contrast'])>=3}
    for j,r in enumerate(records):
        fid=r['exposure_id'];path=files/r['filename'];det=Detector(path,r,c,ee)
        m=pc.metadata(meta,fid.removesuffix('_cal'));raw,x0,y0=pc.wide_images(pixels,m,subtract_template=False)
        cal=fits.getheader(path,'SCI');factor=cal['PHOTMJSR']*cal['PIXAR_SR']*1e12  # microJy/native pixel per DN/s
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning);res=raw-np.nanmedian(raw,axis=0)
        with TemporaryDirectory() as td:
            t=Path(td)
            b=trials(t/'e.csv',map_path,t/'g.npz',np.array([[76418.11928,0,0,0,0,254.746747]]),2460955.765061979)
        take=q['frame']==j;sx,sy=det.native(b['x_path'][0,take],b['y_path'][0,take])
        for i in range(len(raw)):
            if (j,i) not in good:continue
            row=good[j,i];cx=float(row['x']);cy=float(row['y'])
            px=sx[i]-sx[i,2]+cx;py=sy[i]-sy[i,2]+cy
            # Odd local region includes the PSF core, wings and local sky.
            xmin=int(round(cx))-10;ymin=int(round(cy))-10
            cut=np.s_[ymin-y0:ymin-y0+21,xmin-x0:xmin-x0+21]
            fcut=flat[ymin:ymin+21,xmin:xmin+21]
            with np.errstate(invalid='ignore',divide='ignore'):
                img=res[i][cut]*factor/fcut
            kernel=np.zeros((21,21));w=np.array([1,2,2,2,1])/8
            for x,y,wt in zip(np.rint((px-xmin)*8).astype(int),np.rint((py-ymin)*8).astype(int),w,strict=True):
                patch=search.clipped_stamp(kernel.shape,y//8-det.half,x//8-det.half,det.phase_grid[y%8,x%8])
                if patch:
                    ind,psf=patch;kernel[ind]+=wt*psf
            fit=fit_flux(img,kernel)
            if fit is None:continue
            flux,n,support=fit
            yield dict(frame=fid,integration=i+1,mjd_utc=m['times'][i],flux_microjy=flux,good_pixels=n,supported_psf_flux=support)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube1','cube2','map1','map2','fits1','fits2','meta1','pixels1','meta2','pixels2','files','ee','flat','output']:p.add_argument(k,type=Path)
    a=p.parse_args();flat=fits.getdata(a.flat,'SCI');rr=[]
    for c,m,f,h,px in [(a.cube1,a.map1,a.fits1,a.meta1,a.pixels1),(a.cube2,a.map2,a.fits2,a.meta2,a.pixels2)]:
        rr.extend(visit(c,m,f,h,px,a.files,a.ee,flat))
    with a.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rr[0]);writer.writeheader();writer.writerows(rr)
    for prefix in ['jw063','jw094']:
        half=[r for r in rr if r['frame'].startswith(prefix)]
        per=[np.median([s['flux_microjy'] for s in half if s['frame']==v]) for v in sorted({r['frame'] for r in half})]
        print(prefix,len(half), np.round(per,4), 'median',np.median(per),'scaled MAD',1.4826*np.median(abs(per-np.median(per))))


if __name__=='__main__':main()
