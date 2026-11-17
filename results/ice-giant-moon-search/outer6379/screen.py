"""Screen the full native NIRCam detector out to and beyond Oberon's orbit."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage as nd
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
import geometry,ephemeris,prepare


def dump(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def detector(w,at):
    # Interpolate the full SIP WCS on a uniform half-arcsecond tangent grid.
    axis=np.linspace(-52,52,209);e,n=np.meshgrid(axis,axis)
    def exact(e,n):
        return np.array(w.all_world2pix(*prepare.inverse_sky(*at,e/geometry.ARCSEC,n/geometry.ARCSEC),0,quiet=True))
    arrays=exact(e,n)
    def mapping(e,n):
        coords=[(n+52)*2,(e+52)*2]
        return [nd.map_coordinates(table,coords,order=1,mode='constant',cval=-1e5) for table in arrays]
    check=np.linspace(-50,50,81);e,n=np.meshgrid(check+.13,check-.19);e=e.ravel();n=n.ravel()
    err=np.max(np.hypot(*(np.array(mapping(e,n))-exact(e,n))))
    if err>.06:raise ValueError(f'WCS approximation error {err}')
    return mapping,float(err)


def make_frame(file,center,eph):
    with fits.open(file) as h:
        p=h[0].header;s=h['SCI'];w=WCS(s.header,naxis=2);image=s.data.copy();dq=h['DQ'].data.copy()
        start,stop=p['EXPSTART'],p['EXPEND'];times=np.linspace(start,stop,17)
        band=p['FILTER'];visit=p['VISIT_ID']
    _,_,ids,pos=ephemeris.interpolate(eph,times[None,:])
    pos=pos[:,0];b=pos[0]
    ra,dec=w.all_pix2world(float(center['x']),float(center['y']),0)
    native,err=detector(w,(ra,dec))
    good=np.isfinite(image)&((dq&np.uint32(~4 & 0xffffffff))==0)
    yy,xx=np.indices(image.shape,dtype='f4')
    # The entire region interior to the outer rings is excluded in projected distance.
    valid=np.hypot(xx-float(center['x']),yy-float(center['y']))>285
    for k,body in enumerate(ids[1:],1):
        sky=geometry.sky_xy(pos[k,:,0],pos[k,:,1],b[:,0],b[:,1],1.,[0.,0.])
        if np.any(abs(sky.mean(0))>56):continue
        xp,yp=native(-sky[:,0],sky[:,1])
        radius=95 if body in [701,702,703,704,705] else 15 if body in [715,712,714] else 10
        left,top=int(np.floor(np.min(xp)-radius)),int(np.floor(np.min(yp)-radius))
        right,bottom=int(np.ceil(np.max(xp)+radius))+1,int(np.ceil(np.max(yp)+radius))+1
        sl=np.s_[max(0,top):min(image.shape[0],bottom),max(0,left):min(image.shape[1],right)]
        if min(image.shape[0],bottom)<=max(0,top) or min(image.shape[1],right)<=max(0,left):continue
        for x,y in zip(xp,yp):valid[sl]&=(xx[sl]-x)**2+(yy[sl]-y)**2>radius**2
    good&=valid
    def smooth(scale):
        den=nd.gaussian_filter(good.astype('f4'),scale)
        val=nd.gaussian_filter(np.where(good,image,0),scale)/np.maximum(den,.001)
        return val,den
    n,frac=smooth(1.1);back,_=smooth(6.)
    filtered=n-back
    valid &= (frac>.82)&(xx>8)&(yy>8)&(xx<image.shape[1]-9)&(yy<image.shape[0]-9)
    # Local detector tiles calibrate glare, striping and background variation.
    score=np.zeros_like(filtered)
    for y in range(0,image.shape[0],128):
        for x in range(0,image.shape[1],128):
            sl=np.s_[y:y+128,x:x+128];use=valid[sl];v=filtered[sl][use]
            if len(v)<100:valid[sl]=False;continue
            med=np.median(v);scale=1.4826*np.median(abs(v-med))
            if scale<=0:valid[sl]=False;continue
            score[sl]=(filtered[sl]-med)/scale
    middle=len(times)//2
    los,east,north=geometry.tangent_basis(b[middle,0],b[middle,1])
    v1,v2,_=geometry.basis('Uranus',times[middle]+2400000.5)
    proj=np.array([[east@v1,east@v2],[north@v1,north@v2],[los@v1,los@v2]])
    t=Time(times,format='mjd',scale='utc').tdb.jd-b[:,4]*geometry.AU/geometry.C/86400
    return dict(t=t,distance=b[middle,4]*geometry.AU,proj=proj,map=native,score=score,good=valid.astype('f4'),
                wcs_error=err,visit=visit,band=band)


def run(args):
    sources=[(ROOT.parent/'jwst6379_9482/centers_6379.csv',args.long),(ROOT/'short_centers.csv',args.short)]
    # Cubic splines from a public cached ephemeris; each FITS remains native.
    frames=[];files=[]
    for table,eph in sources:
        for r in csv.DictReader(table.open()):
            name=r['exposure_id']+'.fits';file=args.files/name
            f=make_frame(file,r,eph)
            frames.append(f);files.append(dict(filename=name,sha256=hashlib.sha256(file.read_bytes()).hexdigest(),band=f['band'],
                          visit=f['visit'],max_WCS_fit_error_px=round(f['wcs_error'],4)))
    dump('inputs.csv',files)
    epoch=float(np.median([np.median(f['t']) for f in frames]));threshold=1.
    summary=[];hits=[];allmax=[]
    # Phase steps represent up to 1,400 km in the orbital plane. The search is a coarse screen.
    for radius in np.arange(155000.,650001.,2000.):
        nphase=int(np.ceil(2*np.pi*radius/1400))
        phi=np.arange(nphase)*2*np.pi/nphase
        z=[];valid=[]
        for f in frames:
            theta=phi[:,None]+geometry.frequency('Uranus',radius)*(f['t']-epoch)*86400
            u=radius*np.cos(theta);v=radius*np.sin(theta);m=f['proj'];den=f['distance']+m[2,0]*u+m[2,1]*v
            e=(m[0,0]*u+m[0,1]*v)*geometry.ARCSEC/den;n=(m[1,0]*u+m[1,1]*v)*geometry.ARCSEC/den
            x,y=f['map'](e,n)
            val=nd.map_coordinates(f['score'],[y,x],order=1,mode='constant')
            support=nd.map_coordinates(f['good'],[y,x],order=1,mode='constant')
            use=(support>.98).all(1)
            weighted=np.mean(val,axis=1)
            if use.sum()>=60:
                med=np.median(weighted[use]);mad=1.4826*np.median(abs(weighted[use]-med))
                weighted=(weighted-med)/max(.08,mad)
            else:use[:]=False
            z.append(np.where(use,np.clip(weighted,-8,8),0));valid.append(use)
        z,valid=np.array(z),np.array(valid)
        halves=[]
        for sl in [slice(0,5),slice(5,10)]:
            counts=valid[sl].sum(0);halves.append(np.where(counts>=2,z[sl].sum(0)/np.sqrt(np.maximum(1,counts)),-100))
        eligible=(halves[0]>-100)&(halves[1]>-100)&valid[10:].all(0)
        def ranks(values):
            full=values[:10].sum(0)/np.sqrt(np.maximum(1,valid[:10].sum(0)))
            partial=[values[sl].sum(0)/np.sqrt(np.maximum(1,valid[sl].sum(0))) for sl in (slice(0,5),slice(5,10))]
            return np.where(eligible,np.minimum.reduce([full/7,partial[0]/3,partial[1]/3,values[10:].sum(0)/np.sqrt(2)/5,
                                                      values[10:].min(0)/2]),0)
        positive=ranks(z);negative=ranks(-z)
        summary.append(dict(radius_km=int(radius),tracks=nphase,supported=int(eligible.sum()),positive=int((positive>=1).sum()),
                            negative=int((negative>=1).sum()),best=round(float(positive.max()),5)))
        ix=np.flatnonzero(positive>=threshold)
        for k in ix:hits.append(dict(radius_km=radius,phase_deg=np.rad2deg(phi[k]),rank=positive[k],scores=' '.join(f'{v:.2f}' for v in z[:,k])))
        if eligible.any():
            k=positive.argmax();allmax.append(dict(radius_km=radius,phase_deg=np.rad2deg(phi[k]),rank=positive[k]))
    dump('coverage.csv',summary);dump('highest.csv',sorted(allmax,key=lambda r:-r['rank'])[:25])
    if hits:dump('crossings.csv',hits)
    (ROOT/'summary.json').write_text(json.dumps(dict(epoch_emission_jd_tdb=epoch,frames=len(frames),grid=sum(r['tracks'] for r in summary),
           supported=sum(r['supported'] for r in summary),positive=sum(r['positive'] for r in summary),
           inverted=sum(r['negative'] for r in summary)),indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('files','long','short'):ap.add_argument(k,type=Path)
    run(ap.parse_args())
