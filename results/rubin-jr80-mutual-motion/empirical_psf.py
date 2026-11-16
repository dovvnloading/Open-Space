"""Construct per-visit stellar PSFs and fit JR80 using actual reference-star pixels."""
import csv
import gzip
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter, maximum_filter
from scipy.optimize import least_squares
from cutouts import ROOT, DATA, read_hdus, read_cutout, rows

STARS=DATA/"stellar"
def load_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def fit_pair(h, data, variance, psf, radial=7.0):
    yy,xx=np.indices(data.shape); xx=xx-(float(h["CRPIX1"])-1); yy=yy-(float(h["CRPIX2"])-1)
    ok=np.isfinite(data)&np.isfinite(variance)&(variance>0)&(np.hypot(xx,yy)<radial)
    xx,yy,image,sigma=xx[ok],yy[ok],data[ok],np.sqrt(variance[ok])
    sky=np.median(image[np.hypot(xx,yy)>radial-2])
    scale=max(np.sum((image-sky)[np.hypot(xx,yy)<4.5]),1)
    yy0,xx0=(np.array(psf.shape)-1)/2
    source=spline_filter(psf,order=3)
    def profile(x,y):
        return map_coordinates(source,[yy-y+yy0,xx-x+xx0],order=3,mode="constant",prefilter=False)
    def single(v):
        flux,x,y,bg=v
        return (image-bg-scale*flux*profile(x,y))/sigma
    def double(v):
        flux,x,y,bg,dx,dy,p=v
        return (image-bg-scale*flux*((1-p)*profile(x-p*dx,y-p*dy)+p*profile(x+(1-p)*dx,y+(1-p)*dy)))/sigma
    low=[0,-4,-4,-np.inf];high=[10,4,4,np.inf]
    one=least_squares(single,[1.,0.,0.,sky],bounds=(low,high))
    fits=[]
    for a in [0,45,90,135]:
        rad=np.deg2rad(a)
        fits.append(least_squares(double,[*one.x,2.8*np.sin(rad),2.8*np.cos(rad),.49],
                                  bounds=(low+[-6,-6,.1],high+[6,6,.9]),max_nfev=150))
    best=min(fits,key=lambda f:f.cost); v=best.x
    mat=np.array([[float(h[f"PC{i}_{j}"])*float(h[f"CDELT{i}"]) for j in (1,2)] for i in (1,2)])*3600
    xy=mat@v[4:6]
    if xy[0]<0: xy=-xy
    return dict(delta_chi2=max(0,2*(one.cost-best.cost)),east=xy[0],north=xy[1],
                sep=np.linalg.norm(xy),pa=np.degrees(np.arctan2(*xy))%180,
                fraction=min(v[6],1-v[6]))

def aligned_profile(planes):
    (h,a),(_,variance),(_,psf)=planes
    y,x=np.indices(a.shape);xc=float(h["CRPIX1"])-1;yc=float(h["CRPIX2"])-1
    rad=np.hypot(x-xc,y-yc)
    sky=np.nanmedian(a[(rad>9)&(rad<13)]);img=a-sky
    iv=(variance>0)&np.isfinite(img)&np.isfinite(variance)
    if (iv&(rad<8)).sum()<180:
        return None,"pixels"
    # Reject visible companions inside 3–8 pixels with their own local peak.
    peak=np.nanmax(img[rad<3])
    loc=(maximum_filter(np.where(iv,img,-np.inf),size=3)==img)&(img>.18*peak)&(img>5*np.sqrt(variance))
    if np.any(loc&(rad>3)&(rad<8)):
        return None,"neighbour"
    ps=spline_filter(psf); py,px=(np.array(psf.shape)-1)/2; mask=iv&(rad<5.5)
    def res(v):
        amp,cx,cy,bg=v
        pred=amp*map_coordinates(ps,[y[mask]-cy+py,x[mask]-cx+px],order=3,mode="constant",prefilter=False)+bg
        return (img[mask]-pred)/np.sqrt(variance[mask])
    fit=least_squares(res,[max(np.sum(img[(rad<5)&iv]),10),xc,yc,0],bounds=([0,xc-3,yc-3,-np.inf],[np.inf,xc+3,yc+3,np.inf]))
    if np.hypot(fit.x[1]-xc,fit.x[2]-yc)>1.0:
        return None,"offset"
    # Only the centering estimate uses the supplied PSF; shape is measured from
    # star pixels. Shape acceptance uses a circular aperture concentration test.
    new_y,new_x=np.mgrid[-9:10,-9:10]
    blank=np.where(iv,img-fit.x[3],0)
    sample=map_coordinates(blank,[new_y+fit.x[2],new_x+fit.x[1]],order=3,mode="constant")
    rr=np.hypot(new_x,new_y);total=sample[rr<6.5].sum()
    if not (total>0 and .40 < sample[rr<2.5].sum()/total < 1.15):
        return None,"concentration"
    return sample/total,"ok"

def stack(stars):
    # Median of registered normalized pixel grids rejects unshared neighbours.
    model=np.median(np.stack(stars),axis=0)
    yy,xx=np.indices(model.shape);rho=np.hypot(xx-9,yy-9)
    model=model-np.median(model[rho>7])
    model[rho>8.5]=0
    return model/model.sum()

def write(name,data):
    with (ROOT/"output"/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator="\n")
        w.writeheader();w.writerows(data)

def main():
    refs=load_csv(STARS/"candidates.csv")
    available={r["id"]:r for r in load_csv(STARS/"stamps.csv")}
    quality=[];results=[]
    for t in rows("target.csv"):
        tid=t["measurement_id"];slist=[]
        for r in refs:
            if r["target"]!=tid or r["mid"] not in available:continue
            planes=read_hdus(gzip.decompress((STARS/available[r["mid"]]["file"]).read_bytes()))
            pix,flag=aligned_profile(planes)
            quality.append(dict(id=r["mid"],target=tid,status=flag))
            if pix is not None:slist.append((r["mid"],pix,planes))
        if len(slist)<3:
            raise RuntimeError(f"Only {len(slist)} isolated reference stars for {tid}")
        def fit(planes, members):
            (h,img),(_,var),_=planes
            return fit_pair(h,img,var,stack([s[1] for s in members]))
        planes=read_cutout(tid,"Science")
        results.append(dict(id=tid,target=tid,mode="target_all",excluded="",nref=len(slist),**fit(planes,slist)))
        for i,(sid,pix,sp) in enumerate(slist):
            others=slist[:i]+slist[i+1:]
            results.append(dict(id=sid,target=tid,mode="leave_star_out",excluded=sid,nref=len(others),**fit(sp,others)))
            results.append(dict(id=tid,target=tid,mode="target_jackknife",excluded=sid,nref=len(others),**fit(planes,others)))
    write("empirical_star_quality.csv",quality)
    write("empirical_fits.csv",results)
if __name__=="__main__":
    main()
