"""Extract field-star apertures from CHEOPS subarrays using the visit catalogue.

The roll transform follows the CHEOPS PIPE star-catalogue convention
(https://github.com/alphapsa/PIPE/blob/main/pipe/syntstar.py).
"""
from contextlib import contextmanager
from dataclasses import dataclass

import numpy as np
from astropy.io import fits


@dataclass
class Frame:
    cube: object
    idx: object
    good: object
    xc: object
    yc: object
    scale: float
    alignment: dict


def _spread(indices, n):
    return indices[np.linspace(0, len(indices)-1, min(n,len(indices)),dtype=int)]


def _frame(hdus,ref):
    image = next(h for h in hdus if h.header.get("NAXIS")==3 and h.data.ndim==3)
    cube = image.data
    n,height,width = cube.shape
    header = hdus[0].header.copy()
    header.update(image.header)
    mjd = next(h.data["MJD_TIME"] for h in hdus if isinstance(h,fits.BinTableHDU)
               and "MJD_TIME" in h.columns.names and len(h.data)==n)
    rt = np.asarray(ref["MJD_TIME"])
    if len(rt)<2 or not np.all(np.diff(rt)>0):
        raise ValueError("DRP reference cadences are not strictly increasing")
    r = np.clip(np.searchsorted(rt,mjd),1,len(rt)-1)
    idx = np.where(abs(rt[r]-mjd)<abs(rt[r-1]-mjd),r,r-1)
    good = abs(mjd-rt[idx])<0.45*np.median(np.diff(rt))
    if 'EVENT' in ref.names:
        good &= ref['EVENT'][idx]==0
    xc,yc = (np.asarray(ref[k][idx],float) for k in ('CENTROID_X','CENTROID_Y'))
    good &= np.isfinite(xc) & np.isfinite(yc) & np.isfinite(ref['ROLL_ANGLE'][idx])
    if np.count_nonzero(good)<5:
        raise ValueError("Fewer than five valid cadence matches for alignment")
    mx,my = np.median(xc[good]),np.median(yc[good])
    candidates = []
    if 0<=mx<width and 0<=my<height:
        candidates.append(('already local',xc,yc,False))
    if 'X_WINOFF' in header and 'Y_WINOFF' in header:
        candidates.append(('detector minus subarray origin',xc-float(header['X_WINOFF']),
                           yc-float(header['Y_WINOFF']),False))
    if 'CTYPE1' in header and 'RA_TARG' in header and 'DEC_TARG' in header:
        from astropy.wcs import WCS
        try:
            wcs = WCS(header,naxis=2).celestial
            if wcs.has_celestial:
                xw,yw=wcs.world_to_pixel_values(float(header['RA_TARG']),float(header['DEC_TARG']))
                if np.isfinite(xw) and np.isfinite(yw):
                    candidates.append(('celestial WCS target anchor plus measured jitter',xc-mx+xw,yc-my+yw,False))
        except (ValueError,TypeError):
            pass
    if all(k in ref.names for k in ('LOCATION_X','LOCATION_Y','FLUX')):
        # PIPE defines the nominal window centre as offset + half the image size;
        # pycheops uses CENTROID minus LOCATION for the measured target jitter.
        candidates.append(('nominal centre plus DRP centroid minus LOCATION',
                           width/2 + xc-ref['LOCATION_X'][idx],height/2 + yc-ref['LOCATION_Y'][idx],True))
    failures=[]
    for method,x,y,confirm_identity in candidates:
        ok = good & np.isfinite(x) & np.isfinite(y)
        check=_verify_alignment(cube,ref,idx,ok,x,y,confirm_identity)
        if check.get('accepted'):
            check.update(method=method,x_winoff=header.get('X_WINOFF'),y_winoff=header.get('Y_WINOFF'),
                         preceding_rejections=failures)
            return Frame(cube,idx,ok,x,y,float(header.get('PIXSCALE',1.01)),check)
        failures.append((method,check))
    raise ValueError('No image-verified central alignment: '+str(failures))


def _verify_alignment(cube,ref,idx,good,xc,yc,confirm_identity):
    _,height,width=cube.shape
    if np.count_nonzero(good)<5:
        return {'accepted':False,'reason':'insufficient usable coordinates'}
    # Verify the transformed location against a sample of actual target images.
    # Keep the measured floating-point centroid (and its jitter) for every cadence.
    offsets = []
    dominant = []
    flux_ratios = []
    gy,gx = np.mgrid[-24:25,-24:25]
    for k in _spread(np.flatnonzero(good),9):
        ix,iy = round(xc[k]),round(yc[k])
        if not (24<=ix<width-24 and 24<=iy<height-24):
            continue
        stamp = np.asarray(cube[k,iy-24:iy+25,ix-24:ix+25],float)
        sx=gx+ix-xc[k]
        sy=gy+iy-yc[k]
        rr=sx*sx+sy*sy
        live=np.isfinite(stamp) & (stamp!=0)
        sky=stamp[live & (rr>=21**2) & (rr<24**2)]
        aperture=live & (rr<18**2)
        if len(sky)<75 or aperture.sum()<700:
            continue
        bg=np.median(sky)
        noise=max(1.4826*np.median(abs(sky-bg)),np.finfo(float).eps)
        weights=np.where(aperture,stamp-bg,0.)
        net=weights.sum()
        if net>5*noise*np.sqrt(aperture.sum()):
            offsets.append(float(np.hypot((sx*weights).sum(),(sy*weights).sum())/net))
            if confirm_identity:
                # This additional test is confined to a nominal-centre fallback.
                # Require target-like integrated flux and dominance over off-centre peaks.
                from scipy.ndimage import uniform_filter
                im = np.asarray(cube[k],float)
                pixels = np.isfinite(im) & (im!=0)
                baseline = np.nanmedian(im[pixels])
                smoothed = uniform_filter(np.where(pixels,im-baseline,0.),size=15,mode='constant')
                sy_all,sx_all = np.ogrid[:height,:width]
                radius2 = (sx_all-xc[k])**2+(sy_all-yc[k])**2
                centre = float(np.max(smoothed[radius2<15**2]))
                elsewhere = float(np.max(smoothed[radius2>35**2]))
                dominant.append(centre>max(0.,1.25*elsewhere))
                if ref['FLUX'][idx[k]]>0:
                    flux_ratios.append(float(net/ref['FLUX'][idx[k]]))
    if len(offsets)<3:
        return {'accepted':False,'checked_frames':len(offsets),'reason':'weak or out-of-window target image'}
    distance=float(np.median(offsets))
    check=dict(accepted=distance<=6,checked_frames=len(offsets),median_offset_pixels=distance)
    if confirm_identity:
        ratio=float(np.median(flux_ratios)) if flux_ratios else None
        check.update(dominant_frames=sum(map(int,dominant)),median_image_to_drp_flux=ratio)
        check['accepted']=bool(check['accepted'] and len(offsets)>=5 and np.quantile(offsets,.8)<=5
                               and np.mean(dominant)>=.8 and ratio is not None and .25<=ratio<=1.8)
    return check


_Y, _X = np.mgrid[-18:19,-18:19]


def _measure(cube,x,y,indices):
    """Finite aperture pixels and local sky in batches to bound temporary storage."""
    values=[]
    for start in range(0,len(indices),48):
        k=indices[start:start+48]
        ix,iy=np.rint(x[k]).astype(int),np.rint(y[k]).astype(int)
        xx=ix[:,None,None]+_X
        yy=iy[:,None,None]+_Y
        # Coordinates (not the pixel data) keep the actual fractional position.
        dx=xx-x[k,None,None]
        dy=yy-y[k,None,None]
        rr=dx**2+dy**2
        stamp=np.asarray(cube[k[:,None,None],yy,xx],float)
        valid=np.isfinite(stamp) & (stamp!=0)
        ap=rr<100
        sky=valid & (rr>196) & (rr<324)
        ok=(np.sum(sky,axis=(1,2))>=100) & np.all(valid | ~ap,axis=(1,2))
        if not ok.any():
            continue
        k,stamp,sky,ap,dx,dy=(v[ok] for v in (k,stamp,sky,ap,dx,dy))
        sky_vals=np.where(sky,stamp,np.nan).reshape(len(k),-1)
        bg=np.nanmedian(sky_vals,axis=1)
        sigma=1.4826*np.nanmedian(abs(sky_vals-bg[:,None]),axis=1)
        w=np.where(ap,stamp-bg[:,None,None],0.)
        flux=w.sum(axis=(1,2))
        error=np.sqrt(np.maximum(flux,1)+ap.sum(axis=(1,2))*sigma**2)
        delta=np.hypot((w*dx).sum(axis=(1,2)),(w*dy).sum(axis=(1,2)))/np.maximum(flux,1.)
        values.append(np.column_stack((k,flux,error,bg,delta)))
    return np.vstack(values) if values else np.empty((0,5))


def extract(paths, reference, directory, log):
    """Extract resolved catalogue sources showing a measurable core in probe frames."""
    catalogues = [p for p in paths if 'starcatalog' in p.name.lower()]
    images = [p for p in paths if 'subarray' in p.name.lower() and
              any(k in p.name.lower() for k in ('sci_cor','sci_cal'))]
    log.update(catalogues=len(catalogues),image_products=len(images),extracted=0)
    if not (catalogues and images and reference):
        log['status']='missing calibrated subarray, catalogue or DRP reference'
        return
    image=sorted(images,key=lambda p:('sci_cor' not in p.name.lower(),p.name))[0]
    cat=fits.getdata(catalogues[0],1)
    ref=fits.getdata(reference,1)
    refsys=fits.getheader(reference,1).get('TIMESYS',fits.getheader(reference,0).get('TIMESYS','TDB'))
    with fits.open(image,memmap=None) as hdus:
        frame=_frame(hdus,ref)
        cube=frame.cube
        n,height,width=cube.shape
        roll=np.deg2rad(ref['ROLL_ANGLE'][frame.idx])
        co,si=np.cos(roll),np.sin(roll)
        ra,dec=(np.asarray(cat[k],float) for k in ('RA','DEC'))
        dx=((ra[0]-ra+180)%360-180)*np.cos(np.deg2rad(dec[0]))*3600/frame.scale
        dy=(dec-dec[0])*3600/frame.scale
        dist=np.hypot(dx,dy)
        # At least one measurement must keep the full 18-pixel sky ring on-chip.
        x0,y0=np.median(frame.xc[frame.good]),np.median(frame.yc[frame.good])
        furthest=max(np.hypot(x1-x0,y1-y0) for x1 in (18,width-19)
                     for y1 in (18,height-19))+np.max(np.hypot(frame.xc[frame.good]-x0,frame.yc[frame.good]-y0))
        possible=np.isfinite(dist) & (dist>=18) & (dist<=furthest)
        log.update(catalogue_rows=len(cat),near_target=int(np.count_nonzero(dist[1:]<18)),
                   outside=int(np.count_nonzero(~np.isfinite(dist[1:]) | (dist[1:]>furthest))),
                   too_few=0,unresolved_core=0,image=image.name,catalogue=catalogues[0].name,
                   alignment=frame.alignment,geometric_candidates=int(possible.sum()),
                   pilot='up to 15 uniform frames; >=5 usable, median aperture S/N >=3, median core offset <=5 pixels',
                   background='local median 14-18 pixel annulus subtracted from a radius 10 aperture',
                   auxiliaries='roll measured at matched cadence; smear and contamination from primary DRP aperture are proxies')
        for j in np.flatnonzero(possible):
            x=frame.xc+dx[j]*co+dy[j]*si
            y=frame.yc-dx[j]*si+dy[j]*co
            ok=frame.good & (x>=19) & (x<width-19) & (y>=19) & (y<height-19)
            indices=np.flatnonzero(ok)
            if len(indices)<48:
                log['too_few']+=1
                continue
            probe=_measure(cube,x,y,_spread(indices,15))
            if (len(probe)<5 or np.median(probe[:,1]/probe[:,2])<3 or np.median(probe[:,4])>5):
                log['unresolved_core']+=1
                continue
            values=_measure(cube,x,y,indices)
            if len(values)<48 or np.median(values[:,1])<=0:
                log['too_few']+=1
                continue
            ri=frame.idx[values[:,0].astype(int)]
            cols=[fits.Column(name=k,format='D',array=a) for k,a in (
                    ('BJD_TIME',ref['BJD_TIME'][ri]),('FLUX',values[:,1]),('FLUXERR',values[:,2]),
                    ('BACKGROUND',values[:,3]),('ROLL_ANGLE',ref['ROLL_ANGLE'][ri]))]
            for k in ('SMEARING_LC','CONTA_LC','EVENT','STATUS'):
                if k in ref.names:
                    cols.append(fits.Column(name=k,format='D',array=ref[k][ri]))
            hdr=fits.Header()
            hdr['TIMESYS']=refsys
            hdr['OBJECT']='Gaia '+str(cat['ID'][j]).strip()
            hdr['RA_TARG'],hdr['DEC_TARG']=float(ra[j]),float(dec[j])
            path=directory/f'field-{j}.fits'
            fits.HDUList([fits.PrimaryHDU(header=hdr),fits.BinTableHDU.from_columns(cols)]).writeto(path,overwrite=True)
            log['extracted']+=1
            yield path
    log['status']='catalogue apertures extracted within stated spatial and sampling limits'


@contextmanager
def locator(paths,reference):
    """Keep the same visit cube and frame registration open for multiple event checks."""
    if reference is None:
        yield lambda event,coordinate: {"status":"missing_reference"}
        return
    cat_path = next((p for p in paths if "starcatalog" in p.name.lower()),None)
    image = next((p for p in sorted(paths,key=lambda p:"sci_cor" not in p.name.lower())
                  if "subarray" in p.name.lower() and any(s in p.name.lower() for s in ("sci_cal","sci_cor"))),None)
    if not (cat_path and image):
        yield lambda event,coordinate: {"status":"missing_images_or_catalogue"}
        return
    cat = fits.getdata(cat_path,1)
    ref = fits.getdata(reference,1)
    with fits.open(image,memmap=None) as hdus:
        frame = _frame(hdus,ref)
        idx = frame.idx
        t,roll = ref["BJD_TIME"][idx],ref["ROLL_ANGLE"][idx]
        timesys = fits.getheader(reference,1).get("TIMESYS",fits.getheader(reference,0).get("TIMESYS","TDB"))
        if str(timesys).strip().upper()=="TT":
            from astropy.time import Time
            t = Time(t,format="jd",scale="tt").tdb.jd
        yield lambda event,coordinate: _localize(event,coordinate,frame,cat,t,roll,image)


def locate(event, paths, reference, coordinate):
    with locator(paths,reference) as check:
        return check(event,coordinate)


def _localize(event,coordinate,frame,cat,t,roll,image):
    if coordinate is None or None in coordinate:
        return {"status":"missing_source_position"}
    t0,d = event["epoch_bjd"],event["duration_days"]
    cube,ok = frame.cube,frame.good
    out = np.flatnonzero(ok & (abs(t-t0)>d) & (abs(t-t0)<max(1,4*d)))
    if len(out)<5:
        return {"status":"insufficient_reference_frames"}
    dx = ((cat["RA"][0]-coordinate[0]+180)%360-180)*3600*np.cos(np.deg2rad(cat["DEC"][0]))/frame.scale
    dy = (coordinate[1]-cat["DEC"][0])*3600/frame.scale
    co,si = np.cos(np.deg2rad(roll)),np.sin(np.deg2rad(roll))
    xs,ys = frame.xc+dx*co+dy*si, frame.yc-dx*si+dy*co
    ain = np.flatnonzero(ok & (abs(t-t0)<d*0.43))
    diffs,used = [],set()
    for k in ain:
        rd = abs((roll[out]-roll[k]+180)%360-180)
        pos = np.hypot(xs[out]-xs[k],ys[out]-ys[k])
        choices = np.flatnonzero((rd<2) & (pos<.45))
        if not len(choices):
            continue
        m = int(out[choices[np.argmin(rd[choices]+pos[choices])]])
        if m in used:
            continue
        x,y = round(float(xs[k])), round(float(ys[k]))
        if x<21 or y<21 or x+21>=cube.shape[2] or y+21>=cube.shape[1]:
            continue
        stamp = np.asarray(cube[m,y-20:y+21,x-20:x+21],float)-cube[k,y-20:y+21,x-20:x+21]
        if np.isfinite(stamp).all():
            diffs.append(stamp)
            used.add(m)
    if len(diffs)<12:
        return {"status":"insufficient_roll_matched_pairs","pairs":len(diffs)}
    gy,gx = np.mgrid[-20:21,-20:21]
    rr = gx*gx+gy*gy
    measures = []
    for subset in np.array_split(np.array(diffs),2):
        diff = np.median(subset,axis=0)
        sky = diff[rr>=18**2]
        core = (diff-np.median(sky))*(rr<14**2)
        net = np.sum(core)
        noise = 1.4826*np.median(abs(sky-np.median(sky)))*np.sqrt((rr<14**2).sum())
        if net<=0 or noise<=0:
            return {"status":"unstable_pixel_deficit","pairs":len(diffs)}
        measures.append({"snr":float(net/noise),"offset_pixels":float(np.hypot(
            np.sum(core*gx)/net,np.sum(core*gy)/net))})
    good = all(m["snr"]>=6 and m["offset_pixels"]<=2 for m in measures)
    return {"status":"consistent" if good else "unstable_pixel_deficit",
            "pairs":len(diffs),"split_stack":measures,"image":image.name}
