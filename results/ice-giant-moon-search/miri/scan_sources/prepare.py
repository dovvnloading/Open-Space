"""Register JWST calibrated 2-D exposures with SIP WCS on a planet-centered grid."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from geometry import ARCSEC, ring_mask, sky_xy, tangent_basis
from ephemeris import interpolate


def inverse_sky(ra, dec, e, n):
    los, east, north = tangent_basis(ra,dec)
    vec = los + e[...,None]*east + n[...,None]*north
    return np.rad2deg(np.arctan2(vec[...,1],vec[...,0]))%360, np.rad2deg(
        np.arctan2(vec[...,2],np.hypot(vec[...,0],vec[...,1])))


def resample(image, error, dq, x, y, factor, allow_jumps=False):
    """Bilinear science and variance of a linear combination; require all contributors good."""
    h,w = image.shape
    x0,y0 = np.floor(x).astype(int),np.floor(y).astype(int)
    sx,sy = x-x0,y-y0
    valid = np.isfinite(x+y) & (x0>=0)&(y0>=0)&(x0+1<w)&(y0+1<h)
    result,variance = np.zeros(x.shape),np.zeros(x.shape)
    for ix,iy,weight in [(0,0,(1-sx)*(1-sy)),(1,0,sx*(1-sy)),(0,1,(1-sx)*sy),(1,1,sx*sy)]:
        xi,yi = np.clip(x0+ix,0,w-1),np.clip(y0+iy,0,h-1)
        flux,err = image[yi,xi],error[yi,xi]
        dq_good = dq[yi,xi]==0
        if allow_jumps:
            dq_good |= dq[yi,xi]==4
        valid &= (weight<1e-8) | (np.isfinite(flux) & np.isfinite(err) & (err>0) & dq_good)
        result += weight*np.nan_to_num(flux)*factor
        variance += (weight*np.nan_to_num(err)*factor)**2
    return result.astype('f4'),variance.astype('f4'),valid


def prepare(files, ephemeris, output, centers, size, scale, samples, sigma, pad, allow_jumps=False, control_body=None, offset_target=False):
    if size%2 != 1 or samples%2 != 1 or samples<3 or scale<=0 or sigma<=0:
        raise ValueError("positive scale, positive PSF sigma, odd image and odd trail sample count >=3 required")
    center = np.full(2, (size-1)/2)
    yy,xx = np.indices((size,size))
    fit = {}
    sky_fit = {}
    if centers:
        with open(centers) as f:
            for r in csv.DictReader(f):
                fit[r["exposure_id"]] = (float(r["x"]),float(r["y"]))
                if r.get("ra_deg") and r.get("dec_deg"):
                    sky_fit[r["exposure_id"]] = (float(r["ra_deg"]),float(r["dec_deg"]))
    if offset_target and (not fit or set(fit) != set(sky_fit)):
        raise ValueError("offset-target processing requires sky and pixel centers for every exposure")
    times,ids,visits,bands,records = [],[],[],[],[]
    images,variances,valids,footprints = [],[],[],[]
    for file in sorted(map(Path,files)):
        exp = file.stem
        if centers and exp not in fit:
            continue
        with fits.open(file,memmap=False) as hdul:
            obs, sci = hdul[0].header, hdul["SCI"]
            if (obs.get("TELESCOP")!="JWST" or sci.data.ndim!=2 or "VISIT_ID" not in obs
                    or not file.name.endswith("_cal.fits") or sci.header.get("BUNIT")!="MJy/sr"):
                raise ValueError(f"require individual calibrated 2-D JWST MJy/sr exposure: {file}")
            if "A_ORDER" not in sci.header or "B_ORDER" not in sci.header:
                raise ValueError(f"missing SIP distortion model: {file}")
            wcs = WCS(sci.header, naxis=2)
            ref = fit.get(exp)
            if ref is None:
                # This header pair is authoritative only for exposures tracking the planet itself.
                ref = wcs.all_world2pix(sci.header["MT_RA"], sci.header["MT_DEC"],0)
            ra,dec = sky_fit.get(exp, wcs.all_pix2world(*ref, 0))
            if offset_target:
                if not np.isfinite([*ref,ra,dec]).all() or not (0 <= ra < 360 and -90 <= dec <= 90):
                    raise ValueError(f"invalid planet center for {exp}")
                wx,wy=wcs.all_world2pix(ra,dec,0)
                if np.hypot(wx-ref[0],wy-ref[1]) > 0.5:
                    raise ValueError(f"inconsistent sky and pixel planet center for {exp}")
            wr,wd = inverse_sky(ra,dec,-(xx-center[0])*scale/ARCSEC,(yy-center[1])*scale/ARCSEC)
            dx,dy = wcs.all_world2pix(wr,wd,0,quiet=True)
            flux,var,valid = resample(sci.data,hdul["ERR"].data,hdul["DQ"].data,dx,dy,1e6*(scale/ARCSEC)**2,allow_jumps)
            height,width=sci.data.shape
            footprints.append((dx>=0)&(dy>=0)&(dx<width-1)&(dy<height-1))
            images.append(flux);variances.append(var);valids.append(valid)
            ids.append(exp);visits.append(obs["VISIT_ID"])
            bands.append(f'{obs["INSTRUME"]}/{obs["FILTER"]}/{obs.get("PUPIL","")}')
            times.append(np.linspace(obs["EXPSTART"],obs["EXPEND"],samples))
            with file.open("rb") as f:
                digest = hashlib.file_digest(f,"sha256").hexdigest()
            records.append(dict(exposure_id=exp,filename=file.name,sha256=digest,
                                center_native_x=float(ref[0]),center_native_y=float(ref[1]),
                                center_source="external sky+pixel center" if offset_target else ("limb CSV" if centers else "MT_RA/DEC"),
                                allowed_dq=[0,4] if allow_jumps else [0],
                                mjd_start=obs["EXPSTART"],mjd_end=obs["EXPEND"],
                                integration_count=obs.get("NINTS"),effexptm_s=obs.get("EFFEXPTM"),
                                detector=obs["DETECTOR"],target=obs.get("TARGNAME")))
    if not ids or len(set(ids))!=len(ids) or len(set(bands))!=1 or (centers and set(fit)!=set(ids)):
        raise ValueError("require unique exposures in one band, matching the entire centers table")
    times = np.array(times)
    planet, observer, bodies, positions = interpolate(ephemeris,times)
    if observer!="JWST":
        raise ValueError("wrong ephemeris observer")
    valid = np.array(valids)
    mid = samples//2
    body = positions[0]
    known = sky_xy(positions[1:,:,:,0],positions[1:,:,:,1],body[None,:,:,0],body[None,:,:,1],scale,center)
    rings = []
    geometry_valid = np.array(footprints)
    if control_body is not None and control_body not in bodies[1:]:
        raise ValueError("unknown positive control ID")
    for k in range(len(ids)):
        if not offset_target and planet.upper() not in str(records[k]["target"]).upper():
            raise ValueError("tracking target and requested planet differ; inspect pointing")
        ra,dec,diam,_,delta = body[k,mid]
        rings.append(ring_mask(planet,ra,dec,delta,times[k,mid]+2400000.5,scale,center,xx.shape,pad))
        good = geometry_valid[k]
        good &= ~rings[-1] & (np.hypot(xx-center[0],yy-center[1]) > diam/(2*scale)+pad)
        for j,trail in enumerate(known[:,k]):
            if bodies[j+1]==control_body:
                continue
            radius = max(0,positions[j+1,k,:,2].max()/2/scale)+pad
            bx,by = np.min(trail,axis=0)-radius
            tx,ty = np.max(trail,axis=0)+radius
            if tx<0 or ty<0 or bx>=size or by>=size:
                continue
            x0,y0 = np.maximum(0, np.floor([bx,by]).astype(int))
            x1,y1 = np.minimum(size,np.ceil([tx+1,ty+1]).astype(int))
            sl = np.s_[y0:y1,x0:x1]
            for x,y in trail:
                good[sl] &= (xx[sl]-x)**2+(yy[sl]-y)**2>radius**2
        valid[k] &= good
        records[k]["usable_fraction"] = float(np.mean(valid[k]))
        records[k]["unmasked_control_body"] = control_body
    np.savez_compressed(output, image=images,variance=variances,valid=valid,mjd=times[:,mid],
                        exposure_id=ids,visit=visits,psf_sigma_px=sigma,flux_unit="Jy/pixel",
                        band=bands[0],center_xy=center,scale_arcsec=scale,planet=planet,
                        unmasked_control_body=str(control_body),geometry_valid=geometry_valid,known_path=known,
                        known_id=bodies[1:],known_xy=known[:,:,mid,:].transpose(1,0,2),
                        ring_mask=np.array(rings),path_mjd=times,path_ra=body[:,:,0],
                        path_dec=body[:,:,1],path_delta=body[:,:,4],r_au=body[:,mid,3],
                        exptime_s=(times[:,-1]-times[:,0])*86400)
    Path(str(output)+".inputs.json").write_text(json.dumps(records,indent=2)+"\n")


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ephemeris")
    p.add_argument("output")
    p.add_argument("files",nargs="+")
    p.add_argument("--centers")
    p.add_argument("--size",type=int,default=401)
    p.add_argument("--scale",type=float,default=0.031)
    p.add_argument("--samples",type=int,default=65)
    p.add_argument("--psf-sigma",type=float,required=True)
    p.add_argument("--mask-pad",type=float,default=5.)
    p.add_argument("--allow-corrected-jumps",action="store_true")
    p.add_argument("--offset-target",action="store_true",help="use verified planet sky+pixel centers for exposures tracking another target")
    p.add_argument("--control-body",type=int,help="unmask one named moon for positive-control data only")
    a=p.parse_args()
    prepare(a.files,a.ephemeris,a.output,a.centers,a.size,a.scale,a.samples,a.psf_sigma,a.mask_pad,a.allow_corrected_jumps,a.control_body,a.offset_target)
