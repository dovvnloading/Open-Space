"""Extract field-star apertures from CHEOPS subarrays using the visit catalogue.

The roll transform follows the CHEOPS PIPE star-catalogue convention
(https://github.com/alphapsa/PIPE/blob/main/pipe/syntstar.py).
"""
import numpy as np
from astropy.io import fits


def extract(paths, reference, directory, log):
    """Yield one field curve at a time. Record missing inputs and geometric exclusions."""
    catalogues = [p for p in paths if "starcatalog" in p.name.lower()]
    images = [p for p in paths if "subarray" in p.name.lower() and
              any(s in p.name.lower() for s in ("sci_cal", "sci_cor"))]
    log.update(catalogues=len(catalogues), image_products=len(images), extracted=0)
    if not (catalogues and images and reference):
        log["status"] = "missing calibrated subarray, catalogue or DRP reference"
        return
    # Prefer corrected exposures when both stages describe the same frames.
    image = sorted(images, key=lambda p: ("sci_cor" not in p.name.lower(), p.name))[0]
    cat = fits.getdata(catalogues[0], 1)
    ref = fits.getdata(reference, 1)
    with fits.open(image, memmap=False) as hdus:
        cube = next(h.data for h in hdus if h.header.get("NAXIS") == 3 and h.data.ndim == 3)
        n, height, width = cube.shape
        meta = next(h.data for h in hdus if isinstance(h, fits.BinTableHDU)
                    and "MJD_TIME" in h.columns.names and len(h.data) == n)
        rt = np.asarray(ref["MJD_TIME"])
        if len(rt) < 2 or not np.all(np.diff(rt) > 0):
            raise ValueError("Invalid reference cadence")
        t = meta["MJD_TIME"]
        right = np.clip(np.searchsorted(rt, t), 1, len(rt)-1)
        idx = np.where(abs(rt[right]-t) < abs(rt[right-1]-t), right, right-1)
        consistent = abs(rt[idx]-t) < 0.45 * np.median(np.diff(rt))
        roll = np.deg2rad(ref["ROLL_ANGLE"][idx])
        xc, yc = (ref[k][idx] for k in ("CENTROID_X", "CENTROID_Y"))
        # CHEOPS subarray centroid coordinates are local. Reject detector coordinates.
        if not (0 <= np.nanmedian(xc) < width and 0 <= np.nanmedian(yc) < height):
            raise ValueError("DRP centroid outside subarray")
        scale = float(hdus[1].header.get("PIXSCALE", 1.01))
        ra, dec = (np.asarray(cat[k], float) for k in ("RA", "DEC"))
        dx = ((ra[0]-ra+180)%360-180) * np.cos(np.deg2rad(dec[0])) * 3600/scale
        dy = (dec-dec[0])*3600/scale
        log.update(catalogue_rows=len(cat), near_target=0, outside=0, too_few=0,
                   image=image.name, catalogue=catalogues[0].name)
        gy, gx = np.mgrid[-19:20,-19:20]
        for j in range(1, len(cat)):
            distance = np.hypot(dx[j], dy[j])
            if distance < 18:
                log["near_target"] += 1
                continue
            if distance > np.hypot(width, height)/2+20 or not np.isfinite(distance):
                log["outside"] += 1
                continue
            x = xc+dx[j]*np.cos(roll)+dy[j]*np.sin(roll)
            y = yc-dx[j]*np.sin(roll)+dy[j]*np.cos(roll)
            ok = consistent & (x >= 20) & (x < width-20) & (y >= 20) & (y < height-20)
            indices, flux, errors, backgrounds = [], [], [], []
            for k in np.flatnonzero(ok):
                xx, yy = round(float(x[k])), round(float(y[k]))
                cut = np.asarray(cube[k, yy-19:yy+20, xx-19:xx+20],float)
                r2 = (xx+gx-x[k])**2+(yy+gy-y[k])**2
                aperture = cut[r2 < 10**2]
                sky = cut[(r2>14**2) & (r2<18**2)]
                sky = sky[np.isfinite(sky) & (sky!=0)]
                if len(sky)<100 or not np.all(np.isfinite(aperture) & (aperture!=0)):
                    continue
                bg = np.median(sky)
                f = np.sum(aperture-bg)
                sigma = 1.4826*np.median(abs(sky-bg))
                indices.append(k)
                flux.append(f)
                errors.append(np.sqrt(max(f,1)+len(aperture)*sigma**2))
                backgrounds.append(bg)
            if len(indices) < 30 or np.median(flux)<=0:
                log["too_few"] += 1
                continue
            ri = idx[indices]
            cols = [fits.Column(name=k, format="D", array=a) for k,a in
                    (("BJD_TIME",ref["BJD_TIME"][ri]), ("FLUX",flux),
                     ("FLUXERR",errors), ("BACKGROUND",backgrounds), ("ROLL_ANGLE",ref["ROLL_ANGLE"][ri]))]
            for k in ("SMEARING_LC","CONTA_LC","EVENT","STATUS"):
                if k in ref.names:
                    cols.append(fits.Column(name=k,format="D",array=ref[k][ri]))
            table = fits.BinTableHDU.from_columns(cols)
            hdr = fits.Header()
            hdr["OBJECT"] = "Gaia "+str(cat["ID"][j]).strip()
            hdr["RA_TARG"],hdr["DEC_TARG"] = float(ra[j]),float(dec[j])
            hdr["TIMESYS"] = fits.getheader(reference,1).get("TIMESYS",fits.getheader(reference,0).get("TIMESYS","TDB"))
            path = directory/f"field-{j}.fits"
            fits.HDUList([fits.PrimaryHDU(header=hdr),table]).writeto(path,overwrite=True)
            log["extracted"] += 1
            yield path
    log["status"] = "catalogue apertures extracted within stated spatial and sampling limits"


def locate(event, paths, reference, coordinate):
    """Pair in-event subarrays with out-of-event frames at the same roll.

    A centroid is estimated from a background-subtracted stack of differences in
    detector coordinates, recentered on the catalogue source every exposure.
    Independent time blocks must agree. Close companions remain a Gaia check.
    """
    if reference is None or None in coordinate:
        return {"status":"missing_reference"}
    cat_path = next((p for p in paths if "starcatalog" in p.name.lower()),None)
    image = next((p for p in sorted(paths,key=lambda p:"sci_cor" not in p.name.lower())
                  if "subarray" in p.name.lower() and any(s in p.name.lower() for s in ("sci_cal","sci_cor"))),None)
    if not (cat_path and image):
        return {"status":"missing_images_or_catalogue"}
    cat = fits.getdata(cat_path,1)
    ref = fits.getdata(reference,1)
    t0,d = event["epoch_bjd"],event["duration_days"]
    with fits.open(image,memmap=False) as hdus:
        cube = next(h.data for h in hdus if h.header.get("NAXIS")==3 and h.data.ndim==3)
        mjd = next(h.data["MJD_TIME"] for h in hdus if isinstance(h,fits.BinTableHDU)
                   and "MJD_TIME" in h.columns.names and len(h.data)==len(cube))
        rt = ref["MJD_TIME"]
        r = np.clip(np.searchsorted(rt,mjd),1,len(rt)-1)
        idx = np.where(abs(rt[r]-mjd)<abs(rt[r-1]-mjd),r,r-1)
        ok = abs(mjd-rt[idx])<0.45*np.median(np.diff(rt))
        if "EVENT" in ref.names:
            ok &= ref["EVENT"][idx]==0
        t,roll = ref["BJD_TIME"][idx],ref["ROLL_ANGLE"][idx]
        timesys = fits.getheader(reference,1).get("TIMESYS",fits.getheader(reference,0).get("TIMESYS","TDB"))
        if str(timesys).strip().upper()=="TT":
            from astropy.time import Time
            t = Time(t,format="jd",scale="tt").tdb.jd
        out = np.flatnonzero(ok & (abs(t-t0)>d) & (abs(t-t0)<max(1,4*d)))
        if len(out)<5:
            return {"status":"insufficient_reference_frames"}
        dx = ((cat["RA"][0]-coordinate[0]+180)%360-180)*3600*np.cos(np.deg2rad(coordinate[1]))/1.01
        dy = (coordinate[1]-cat["DEC"][0])*3600/1.01
        co,si = np.cos(np.deg2rad(roll)),np.sin(np.deg2rad(roll))
        xs,ys = ref["CENTROID_X"][idx]+dx*co+dy*si, ref["CENTROID_Y"][idx]-dx*si+dy*co
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
