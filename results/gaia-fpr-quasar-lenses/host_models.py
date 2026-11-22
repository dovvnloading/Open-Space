"""Compare central extended light with a pair of active nuclei and two hosts.

All hypotheses have the same two Gaia point positions, common WCS shift, fitted
PSF, sloping sky, pixel mask and weights. The paired-host alternative gives
each nucleus a separate nonnegative extended flux, scale, axis ratio and PA.
The scores describe the available pixels; they do not measure a redshift.
"""
import hashlib
import io
import json
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image, ImageDraw
from scipy.optimize import least_squares
from scipy.signal import fftconvolve

from imaging import retrieve, moffat


def original_band(evidence, band_key, band_info, received):
    survey, band = band_key.split("/")
    width = round(band_info["psf_search_width_arcsec"]/band_info["pixel_scale_arcsec"])
    for row in evidence["imaging_inputs"]:
        url = row["url"]
        params = parse_qs(urlparse(url).query)
        if str(width) not in params.get("size", []):
            continue
        if survey=="PS1":
            name = params.get("red", [""])[0]
            if not (urlparse(url).hostname=="ps1images.stsci.edu" and f".stk.{band}." in name):
                continue
        elif survey.startswith("Legacy-"):
            if "legacysurvey.org" not in (urlparse(url).hostname or "") or band not in params.get("bands",[""])[0]:
                continue
        else:
            continue
        payload = retrieve(url, {}, received)
        receipt = received[-1]
        receipt["reference_sha256"] = row["sha256"]
        receipt["same_fits_bytes"] = receipt["sha256"]==row["sha256"]
        with fits.open(io.BytesIO(payload), memmap=False) as hdus:
            hdu = hdus[0]
            data = np.asarray(hdu.data,dtype=float)
            if data.ndim==3:
                bands=str(hdu.header.get("BANDS",params["bands"][0])).strip().replace(" ","")
                data=data[bands.index(band)]
            wcs = WCS(hdu.header,naxis=2).celestial
        if data.ndim!=2 or abs(data.shape[1]-width)>1:
            raise ValueError("unexpected FITS band/size")
        return data,wcs
    raise ValueError("no matching primary FITS receipt")


def compare_bands(pixels, image, info):
    scale=info["pixel_scale_arcsec"]
    sep=np.linalg.norm(pixels[0]-pixels[1])
    h=int(np.ceil(sep/2+3.3/scale))
    center=np.rint(pixels.mean(axis=0)).astype(int)
    cx,cy=center
    positions=pixels-center
    sky,noise=info["background_sky_per_pixel"],info["sky_rms_per_pixel"]
    yy,xx=np.indices(image.shape)
    blank=image[np.hypot(xx-pixels[:,0].mean(),yy-pixels[:,1].mean())>h+3]
    if abs(np.nanmedian(blank)-sky)>.05*noise:
        raise ValueError("retrieved sky fails comparison with deposited stamp")
    new_noise=1.4826*np.nanmedian(abs(blank-np.nanmedian(blank)))
    if abs(new_noise/noise-1)>.05:
        raise ValueError("retrieved variance fails comparison with deposited stamp")
    raw=(image[cy-h:cy+h+1,cx-h:cx+h+1]-sky)/noise
    y,x=np.mgrid[-h:h+1,-h:h+1]
    use=np.isfinite(raw)&(x*x+y*y<=h*h)
    if use.sum()<.90*np.pi*h*h:
        raise ValueError("missing central pixels")
    err=np.hypot(1,.035*np.maximum(np.nan_to_num(raw),0))
    corr=info["variance_inflation"]
    f0=info["empirical_fwhm_arcsec"]/scale
    if not f0>0:
        raise ValueError("no empirical field PSF")
    # Pad before convolution so a broad off-centre component can scatter light
    # into the fitting aperture from outside of it.
    pad=max(10,int(3*f0))
    gy,gx=np.mgrid[-h-pad:h+pad+1,-h-pad:h+pad+1]
    ky,kx=np.mgrid[-pad:pad+1,-pad:pad+1]
    def radius2(xx,yy,q,phi):
        c,s=np.cos(phi),np.sin(phi)
        return (c*xx+s*yy)**2+(c*yy-s*xx)**2/q**2
    def forward(params, family, profile):
        a,b,dx,dy,width,pq,pphi,bg,sx,sy=params[:10]
        # Area-equivalent size keeps the integral of the Moffat PSF at unity.
        psf=lambda xx,yy: moffat(pq*radius2(xx,yy,pq,pphi),width)
        out=bg+sx*x+sy*y
        for flux,(px,py) in zip((a,b),positions):
            out+=flux*psf(x-px-dx,y-py-dy)
        kern=psf(kx,ky)
        kern/=kern.sum()
        host=np.zeros_like(gx,dtype=float)
        groups=params[10:]
        if family=="central":
            centers=[groups[:2]]
            profiles=[groups[2:]]
        else:
            centers=positions+[dx,dy]
            profiles=groups.reshape(2,4)
        for (hx,hy),(flux,rad,q,angle) in zip(centers,profiles):
            # Nine subpixel samples suppress cusp/pixel-phase sensitivity.
            for oy in (-1/3,0.,1/3):
                for ox in (-1/3,0.,1/3):
                    t=radius2(gx+ox-hx,gy+oy-hy,q,angle)/rad**2
                    light=np.exp(-.5*t) if profile=="Gaussian" else np.exp(-np.sqrt(t))
                    host+=flux*light/(18*np.pi*q*rad**2)
        out+=fftconvolve(host,kern,mode="same")[pad:-pad,pad:-pad]
        return out

    cap=max(100,raw[use].max()*np.pi*(4*f0)**2)
    pflux=np.asarray(info["point_flux_in_native_units"])/noise
    lo=[0,0,-.6/scale,-.6/scale,.85*f0,.7,-np.pi,-5,-.4,-.4]
    hi=[4*cap,4*cap,.6/scale,.6/scale,1.15*f0,1,np.pi,5,.4,.4]
    initial=[*np.clip(pflux,1,3*cap),0,0,f0,.93,0.,0,0,0]
    hostlo=[0,.20/scale,.35,-np.pi]
    hosthi=[12*cap,3.2/scale,1.,np.pi]
    panel, result={},{}
    for profile in ("Gaussian","exponential"):
        packed={}
        for family in ("central","paired_hosts"):
            g0=positions.mean(axis=0)
            window=.45*sep
            lower=lo+([*(g0-window),*hostlo] if family=="central" else hostlo*2)
            upper=hi+([*(g0+window),*hosthi] if family=="central" else hosthi*2)
            def residual(p):
                return ((forward(p,family,profile)-raw)/err)[use]
            best=None
            for theta, width_arcsec in ((-.8,.6),(.25,1.2),(1.2,.85)):
                scale_guess=width_arcsec/scale
                host_start=[min(cap,float(sum(pflux))*.2),scale_guess,.75,theta]
                if family=="central":
                    tail=[*g0,*host_start]
                else:
                    tail=host_start+[host_start[0],scale_guess*.7,.82,-theta]
                fit=least_squares(residual,[*initial,*tail],bounds=(lower,upper),x_scale="jac",
                                  max_nfev=120,ftol=3e-6,xtol=3e-6)
                if best is None or fit.cost<best.cost:
                    best=fit
            p=best.x
            model=forward(p,family,profile)
            chi=2*best.cost
            n_eff=use.sum()/corr
            bic=chi/corr+len(p)*np.log(n_eff)
            hostshape=p[12:].reshape(1,4) if family=="central" else p[10:].reshape(2,4)
            packed[family]=dict(score=float(bic),chi2_raw=float(chi),parameters=len(p),
                max_weighted_residual=float(np.max(abs(residual(p)))),
                optimizer_converged=bool(best.success),point_flux_native=(p[:2]*noise).tolist(),
                common_shift_arcsec=(p[2:4]*scale).tolist(),fitted_psf_fwhm_arcsec=float(p[4]*scale),
                fitted_psf_q=float(p[5]),host_flux_native=(hostshape[:,0]*noise).tolist(),
                host_major_scale_arcsec=(hostshape[:,1]*scale).tolist(),host_axis_ratios=hostshape[:,2].tolist(),
                fitted_light_center_rel_pixel=(p[10:12].tolist() if family=="central" else None))
            panel[profile+" "+family]=(raw-model)/err
        result[profile]=dict(fits=packed,host_pair_minus_central_score=float(
            packed["paired_hosts"]["score"]-packed["central"]["score"]))
    return dict(pixels_fitted=int(use.sum()),effective_pixels=float(use.sum()/corr),
                covariance_factor=corr,profile_comparisons=result), panel


def montage(panels,path):
    im=Image.new("L",(220*2,145*2),24)
    draw=ImageDraw.Draw(im)
    for i,(label,data) in enumerate(panels.items()):
        x,y=(i%2)*220,(i//2)*145
        draw.text((x+3,y+2),label,fill=240)
        # Fixed mapping for comparing signed residuals from the competing fits.
        a=np.clip(127+35*np.nan_to_num(data,nan=0),0,255)[::-1].astype("uint8")
        im.paste(Image.fromarray(a).resize((126,126)),(x+45,y+16))
    im.save(path)


def run(evidence,out):
    evidence,out=Path(evidence),Path(out)
    out.mkdir(parents=True,exist_ok=True)
    packet_path=evidence/"review_packets.jsonl"
    by_source={}
    for l in packet_path.read_text().splitlines():
        p=json.loads(l)
        if p["prior"] or not p["imaging"]["imaging_model_pass"] or len(p["components"])!=2:
            continue
        ids=tuple(sorted(str(c["source_id"]) for c in p["components"]))
        if ids not in by_source:
            by_source[ids]=p
    systems=[]
    for packet in by_source.values():
        ident=packet["system_id"]
        image_record=packet["imaging"]
        components=packet["components"]
        coord=SkyCoord([c["ra"] for c in components],[c["dec"] for c in components],unit="deg")
        record=dict(id=ident,bands={},inputs=[],errors={})
        for band_key in image_record["matched_bands"]:
            try:
                info=image_record["bands"][band_key]
                data,wcs=original_band(image_record,band_key,info,record["inputs"])
                world=np.column_stack(wcs.world_to_pixel(coord))
                fit,panels=compare_bands(world,data,info)
                name=ident+"-"+band_key.replace("/","-")+"-residuals.png"
                montage(panels,out/name)
                fit["residual_montage"]=name
                record["bands"][band_key]=fit
            except Exception as exc:
                record["errors"][band_key]=f"{type(exc).__name__}: {exc}"
        record["both_filters_tested"]=len(record["bands"])==2
        systems.append(record)
    (out/"host_comparisons.json").write_text(json.dumps(dict(systems=systems,
        review_packet_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
        method="Two PSFs with common shifts, free elliptical PSF and sloping sky in both hypotheses. "
        "Central: one free-centroid elliptical extended component. Paired hosts: independent nonnegative "
        "flux, scale, ellipticity and orientation at each point image. Repeat with Gaussian and exponential profiles. "
        "Positive paired-host-minus-central score favours central light. Scores use the previously measured blank-sky "
        "variance inflation and effective-pixel log penalty, not calibrated posterior odds. Gaussian size is major-axis "
        "sigma; exponential size is major-axis scale length. Residual montages use an identical scale.",
        limits="Seeing, PSF mismatch, bright knots, asymmetric coalescing hosts and projected galaxies remain ambiguous; "
        "a preferred central component provides no foreground redshift and no evidence of a repeated quasar spectrum.",
        ),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(*sys.argv[1:])
