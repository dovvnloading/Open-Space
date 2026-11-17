"""Apply the native RATEINTS confirmation rule, with no implanted flux, to every retained link."""
import argparse,csv
from pathlib import Path
import sys,json
from tempfile import TemporaryDirectory
import numpy as np
from scipy import ndimage
import warnings
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'jwst9482'));sys.path.insert(0,str(base))
from link_visits import geometry,search,trials
from prepare import inverse_sky
import positive_control as pc


def measure(cube,mp,meta,pixels,elem,epoch,tmp):
    c=search.load_npz(cube);nc=search.load_npz(mp)
    b=trials(tmp/'e.csv',mp,tmp/'b.npz',elem,epoch)
    records=json.loads(Path(str(cube)+'.inputs.json').read_text())
    x=b['x'];y=b['y']
    yy,xx=np.mgrid[-20:21,-20:21]
    ann=(xx**2+yy**2>=8**2)&(xx**2+yy**2<=20**2);core=xx**2+yy**2<=3**2
    use=ann|core;ox=xx[use];oy=yy[use];inann=ann[use];incore=core[use]
    # Use native detector-coordinate translation across a visit, as in planted verification.
    data=[]
    for j,r in enumerate(records):
        m=pc.metadata(meta,r['exposure_id'].removesuffix('_cal'))
        im,x0,y0=pc.wide_images(pixels,m)
        ra,dec=m['wcs'].all_pix2world(r['center_native_x'],r['center_native_y'],0)
        s=nc['frame']==j;sc=c['scale_arcsec']/geometry.ARCSEC
        wr,wd=inverse_sky(ra,dec,-(x[:,s]-c['center_xy'][0])*sc,(y[:,s]-c['center_xy'][1])*sc)
        dx,dy=m['wcs'].all_world2pix(wr,wd,0,quiet=True)
        for i,frame in enumerate(im):
            ok=np.isfinite(frame).astype(float);raw=np.nan_to_num(frame)
            def smooth(sigma):return ndimage.gaussian_filter(raw,sigma)/np.maximum(ndimage.gaussian_filter(ok,sigma),.01)
            filtered=smooth(.9)-smooth(4)
            coords=[dy[:,i,None]-y0+oy,dx[:,i,None]-x0+ox]
            val=ndimage.map_coordinates(filtered,coords,order=1,mode='constant',cval=np.nan)
            valid=ndimage.map_coordinates(ok,coords,order=1,mode='constant')>=.9
            data.append(np.where(valid,val,np.nan).astype('f4'))
    arr=np.array(data);mid=len(arr)//2
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        images=[]
        for block in [arr,arr[:mid],arr[mid:]]:
            z=np.nanmedian(block,axis=0);bg=np.nanmedian(z[:,inann],axis=1)
            mad=1.4826*np.nanmedian(abs(z[:,inann]-bg[:,None]),axis=1)
            images.append((z-bg[:,None])/mad[:,None])
    res=[]
    for i,im in enumerate(images[0]):
        k=np.argmax(np.where(incore&np.isfinite(im),im,-np.inf))
        res.append([im[k],ox[k],oy[k],images[1][i,k],images[2][i,k]])
    return res


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube1','cube2','map1','map2','meta1','pixels1','meta2','pixels2','links','veto','output']:p.add_argument(k,type=Path)
    a=p.parse_args()
    rr=list(csv.DictReader((a.links/'passing_pairs.csv').open()))
    keep=[int(r['input_row']) for r in csv.DictReader((a.veto/'link_distances.csv').open())
          if r['early_sign']==r['late_sign']=='1' and r['veto']=='False']
    elem=np.zeros((len(keep),6));elem[:,0]=[float(rr[k]['a_km']) for k in keep];elem[:,-1]=[float(rr[k]['phase_first_deg']) for k in keep]
    epoch=json.loads((a.links/'parameters.json').read_text())['epoch1_tdb_jd']
    with TemporaryDirectory() as td:
        t=Path(td)
        stats=[measure(c,m,h,px,elem,epoch,t) for c,m,h,px in
               [(a.cube1,a.map1,a.meta1,a.pixels1),(a.cube2,a.map2,a.meta2,a.pixels2)]]
    with a.output.open('w',newline='') as f:
        w=csv.writer(f);fields=['peak','dx','dy','first_half','second_half']
        w.writerow(['row','passed',*['feb_'+k for k in fields],*['oct_'+k for k in fields]])
        for i,k in enumerate(keep):
            z=np.array([stats[j][i] for j in range(2)])
            passed=(z[:,0]>=5).all() and (z[:,3:]>=2).all()
            w.writerow([k,passed,*z.ravel()])


if __name__=='__main__':main()
