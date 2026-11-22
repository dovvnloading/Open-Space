"""Forced aperture photometry in the calibrated native Legacy Surveys DR9/DR10 coadds.

The science images are in AB nanomaggies/pixel. Their inverse-variance and mask
images share the grid. Local stellar curves of growth calibrate a smoothing
kernel to the cube seeing model. Both the continuum and the line are aperture
fluxes. No source-independent correction to a total galaxy flux is available.
Data model: https://www.legacysurvey.org/dr10/files/#coadd and its DR9 counterpart.
"""
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from pathlib import Path
import math

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import gaussian_filter

from legacy_bricks import _fetch, _paths


def _hdu(hdul):
    for h in hdul:
        if h.header.get("NAXIS") == 2 and WCS(h.header).celestial.pixel_n_dim == 2:
            return h
    raise ValueError("missing two-dimensional celestial image")


def _cut(h, x, y, size, clip=False):
    ix, iy = int(math.floor(x+.5)), int(math.floor(y+.5))
    x0, y0 = ix-size, iy-size
    if not (0 <= ix < h.shape[1] and 0 <= iy < h.shape[0]):
        raise ValueError("position off image")
    if not clip and (x0 < 0 or y0 < 0 or ix+size >= h.shape[1] or iy+size >= h.shape[0]):
        raise ValueError("aperture or sky box at coadd edge")
    x0, y0 = max(0,x0), max(0,y0)
    x1, y1 = min(h.shape[1],ix+size+1), min(h.shape[0],iy+size+1)
    return np.asarray(h.section[y0:y1, x0:x1]), x-x0, y-y0


def _circle(size, dx, dy, radius):
    """Fractional coverage using 5x5 subdivisions; native pixels are ~0.26 arcsec."""
    yy, xx = np.mgrid[-size:size+1, -size:size+1]
    weights = np.zeros(xx.shape)
    shifts = (np.arange(5)+.5)/5-.5
    for sy in shifts:
        for sx in shifts:
            weights += (xx+sx-dx)**2+(yy+sy-dy)**2 <= radius**2
    return weights / 25.


def _catalog(file, wcs, band):
    with fits.open(file, memmap=True) as hd:
        t = hd[1].data
        cols = {k.lower(): k for k in t.names}
        x, y = wcs.world_to_pixel_values(t["ra"], t["dec"])
        radius = np.array(t[cols["shape_r"]]) if "shape_r" in cols else np.zeros(len(t))
        radius = np.maximum(1.7, 3*np.nan_to_num(radius))
        fs = np.array(t["flux_"+band])
        iv = np.asarray(t["flux_ivar_"+band])
        skyobj = np.asarray(t["type"]).astype(str)
        exists = np.isfinite(x+y) & (np.char.strip(skyobj) != "DUP")
        clean = (exists & (np.char.strip(skyobj) == "PSF") & np.asarray(t["brick_primary"], bool)
                 & (t["maskbits"] == 0) & (t["nobs_"+band] > 0) & (t["fracflux_"+band] < .03)
                 & (t["fracmasked_"+band] < .03) & (fs*np.sqrt(np.maximum(0,iv)) > 80))
        return x[exists], y[exists], radius[exists], clean[exists], fs[exists]


def _curve_kernel(h, cat, cx, cy, scale, seeing, area_seeing, beta):
    """Normalize stellar flux within 6 arcsec and fit three enclosed-light fractions."""
    xx, yy, sizes, good, _ = cat
    dist = np.hypot(xx-cx, yy-cy)*scale
    indices = np.flatnonzero(good & (dist < 440) & (dist > 9))
    indices = indices[np.argsort(dist[indices])]
    sigmas = np.linspace(0, max(1.2, area_seeing/2), 16)
    rgrid = (0.9, 1.5, 2.5)
    alpha = area_seeing/(2*math.sqrt(2**(1/beta)-1))
    ee = lambda radius: 1-(1+(radius/alpha)**2)**(1-beta)
    target = np.array([ee(r)/ee(6) for r in rgrid])
    curves = []
    radius = math.ceil((8+4*max(sigmas))/scale)
    for i in indices[:180]:
        dstar = np.hypot(xx-xx[i], yy-yy[i])*scale
        if np.any((dstar > .01) & (dstar < 7+sizes)):
            continue
        try:
            image, sx, sy = _cut(h["image"],xx[i],yy[i],radius)
            iv, _, _ = _cut(h["invvar"],xx[i],yy[i],radius)
            mask, _, _ = _cut(h["maskbits"],xx[i],yy[i],radius)
        except ValueError:
            continue
        y, x = np.indices(image.shape)
        rr = np.hypot(x-sx,y-sy)*scale
        valid = np.isfinite(image) & np.isfinite(iv) & (iv>0) & (mask==0)
        sky = valid & (rr>6.5) & (rr<8)
        if np.any(~valid & (rr<6+3*max(sigmas))) or sky.sum()<50:
            continue
        data = np.where(valid,image-np.median(image[sky]),0.)
        masks = [_circle(radius,sx-radius,sy-radius,r/scale) for r in (*rgrid,6.)]
        grid = []
        for sig in sigmas:
            smooth = gaussian_filter(data, sig/scale, mode="constant", truncate=3.5) if sig else data
            sums = np.array([np.sum(m*smooth) for m in masks])
            grid.append(sums[:-1]/sums[-1] if sums[-1] > 0 else [np.nan]*3)
        if np.isfinite(grid).all():
            curves.append(grid)
        if len(curves)>=10:
            break
    fwhm_kernel = math.sqrt(max(area_seeing**2-seeing**2,0))/2.35482
    out = dict(method="three local stellar encircled fractions vs cube Moffat", stars=len(curves),
               radius_arcsec=list(rgrid), normalization_radius_arcsec=6, cube_moffat_beta=beta,
               cube_seeing_arcsec=area_seeing, target_fraction=target.tolist(), psf_matched=False)
    if len(curves)<3:
        return fwhm_kernel, dict(out, method="Gaussian FWHM approximation; insufficient clean stars")
    curves = np.asarray(curves)
    med = np.median(curves,axis=0)
    k = int(np.argmin(np.sum((med-target)**2,axis=1)))
    scatter = 1.4826*np.median(np.abs(curves[:,k,:]-med[k]),axis=0)
    # Reproducible core-and-wing tolerances, absolute enclosed-light fractions.
    out.update(observed_fractions=med[k].tolist(), stellar_scatter=scatter.tolist(),
               aperture_fractional_uncertainty=float(max(.08,(abs(med[k,1]-target[1])+scatter[1])/target[1]+.05)),
               psf_matched=bool(np.max(np.abs(med[k]-target))<=.065 and scatter.max()<=.07))
    return float(sigmas[k]),out


def _measure(images, x, y, xcat, ycat, radii, scale, radius, sigma):
    image, ivar, mask, nexp = images
    half = int(math.ceil((radius+3.5*sigma)/scale)+2)
    weights = _circle(half,x-int(math.floor(x+.5)),y-int(math.floor(y+.5)),radius/scale)
    if sigma>0:
        weights = gaussian_filter(weights,sigma/scale,mode="constant",truncate=3.5)
    use = weights>1e-5
    weights = weights[use]
    def aperture(ix,iy):
        z = np.s_[iy-half:iy+half+1,ix-half:ix+half+1]
        if ix-half<0 or iy-half<0 or ix+half>=ivar.shape[1] or iy+half>=ivar.shape[0]:
            return None
        f,v,m,n = (a[z][use] for a in images)
        if not (np.isfinite(v).all() and np.all(v>0) and np.isfinite(f).all() and np.all(m==0) and np.all(n>0)):
            return None
        return float(np.dot(f,weights)), float(np.sqrt(np.sum(weights**2/v)))
    ix, iy = int(math.floor(x+.5)),int(math.floor(y+.5))
    source = aperture(ix,iy)
    if source is None:
        return dict(status="masked_or_uncovered",eligible=False)
    # Fixed nonoverlapping positions; the grids and source veto are independent of target flux.
    step = int(math.ceil(2*(radius+3.5*sigma+.4)/scale))
    blank = []
    offsets = []
    for dy in range(-4,5):
        for dx in range(-4,5):
            ox, oy = dx*step,dy*step
            if not (max(8.,2*(radius+3.5*sigma)) <= math.hypot(ox,oy)*scale <= 55):
                continue
            if np.any(np.hypot(xcat-x-ox,ycat-y-oy)*scale < radii+radius+2*sigma):
                continue
            ap = aperture(ix+ox,iy+oy)
            if ap is not None and .60<=ap[1]/source[1]<=1.65:
                blank.append(ap)
                offsets.append([ox,oy])
    if len(blank)<16:
        return dict(status="insufficient matched blank apertures",eligible=False,blank_apertures=len(blank))
    b = np.asarray(blank)
    sky = np.median(b[:,0])
    residual = (b[:,0]-sky)/b[:,1]
    sd = max(1.,1.4826*np.median(abs(residual)), (np.quantile(residual,.84)-np.quantile(residual,.16))/2)
    sky_error = 1.2533*np.median(b[:,1])*sd/math.sqrt(len(blank))
    # Retain the full local sky correction as an extra systematic, important near faint diffuse light.
    err = float(np.sqrt((source[1]*sd)**2+sky_error**2+sky**2))
    return dict(status="measured",flux=float(source[0]-sky),flux_err=err,error=err,
                raw_flux_nanomaggy=source[0],sky_aperture_nanomaggy=float(sky),blank_noise_scale=float(sd),
                source_formal_error=source[1],blank_apertures=len(blank),blank_offsets_pixels=offsets,
                blank_flux_error_nanomaggy=b.tolist(), aperture_radius_arcsec=radius,
                mask_rule="all optical maskbits zero; positive NEXP and inverse variance at every weighted native pixel")


def _band(rec, cone, cache_dir, cube, band):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    root = Path(getattr(cache_dir,"root",cache_dir)) / "legacy-tractor"
    choices = [p for p in cone.get("coverage",[]) if p.get("contains_center") and p.get("status")=="queried"
               and p["release"]==cone.get("selected_release") and p["cap"]==cone.get("selected_cap")]
    fail = dict(status="unavailable",eligible=False,band=band,receipts=[])
    if not choices:
        return dict(fail,reason="no verified containing native brick")
    entry = choices[0]
    rel, cap, brick = entry["release"], entry["cap"],entry["brickname"]
    tractor,base,_ = _paths(rel,cap,brick)
    names = dict(tractor=tractor)
    for kind in ("image","invvar","psfsize","maskbits","nexp"):
        suffix = kind if kind=="maskbits" else kind+"-"+band
        names[kind] = base+f"legacysurvey-{brick}-{suffix}.fits.fz"
    with ThreadPoolExecutor(max_workers=4) as pool:
        products = dict(zip(names,pool.map(lambda path:_fetch(root,path,192<<20),names.values())))
    provenance = dict(band=band,release=rel,cap=cap,brickname=brick,receipt_keys=names,
                      receipts=[r for _,r in products.values()])
    if not all(p for p,_ in products.values()):
        return {**fail,**provenance,"reason":"missing native coadd products"}
    try:
        ra,dec = float(rec["ra"]),float(rec["dec"])
        with ExitStack() as stack:
            h = {k:_hdu(stack.enter_context(fits.open(p,memmap=False))) for k,(p,_) in products.items() if k!="tractor"}
            unit = str(h["image"].header.get("BUNIT","")).lower().strip()
            if unit and "nanomagg" not in unit and "nmaggy" not in unit and "nmgy" not in unit:
                raise ValueError("inconsistent coadd BUNIT "+unit)
            w = WCS(h["image"].header).celestial
            cx,cy = w.world_to_pixel_values(ra,dec)
            scale = float(np.mean(proj_plane_pixel_scales(w))*3600)
            if not .20<scale<.35:
                raise ValueError("unexpected native pixel scale")
            for tile in h.values():
                wx = WCS(tile.header).celestial
                xx,yy = wx.world_to_pixel_values(ra,dec)
                if tile.shape!=h["image"].shape or math.hypot(cx-xx,cy-yy)>.05:
                    raise ValueError("native coadd grids differ")
            fw,_,_= _cut(h["psfsize"],cx,cy,1)
            seeing = float(np.median(fw))
            if not .35<seeing<4:
                raise ValueError("invalid measured local PSF FWHM")
            points=_catalog(products["tractor"][0],w,band)
            sx,sy,radii,_,_=points
            target_seeing = (cube or {}).get("seeing_arcsec")
            beta = float((cube or {}).get("two_psf_beta") or 2.5)
            kernel, psf = (0.,dict(psf_matched=False,method="no independent cube seeing",stars=0))
            if target_seeing and .5<float(target_seeing)<5 and 1<beta<8:
                kernel,psf=_curve_kernel(h,points,cx,cy,scale,seeing,float(target_seeing),beta)
            half = math.ceil(67/scale)
            cut = [_cut(h[k],cx,cy,half,clip=True) for k in ("image","invvar","maskbits","nexp")]
            arrays = [s[0] for s in cut]
            x,y = cut[0][1:]
            act = np.hypot(sx-cx,sy-cy)*scale
            neighbours = ((act>1.25) & (act<radii+1.5+2*kernel))
            ncenter = int(np.sum(act<=1.25))
            isolated = not np.any(neighbours) and ncenter<=1
            result = _measure(arrays,x,y,sx-cx+x,sy-cy+y,radii,scale,1.5,kernel)
            broad = _measure(arrays,x,y,sx-cx+x,sy-cy+y,radii,scale,3.,kernel)
            cube_sep = None
            if cube and cube.get("ra") is not None and cube.get("dec") is not None:
                cube_sep=float(SkyCoord(ra*u.deg,dec*u.deg).separation(
                    SkyCoord(float(cube["ra"])*u.deg,float(cube["dec"])*u.deg)).arcsec)
            result.update(provenance)
            result.update(unit="nanomaggies",native_bunit=unit or "unspecified; official coadd data model",
                          ra=ra,dec=dec,pixel_scale_arcsec=scale,seeing_arcsec=seeing,
                          n_central_counterparts=ncenter, nearby_contaminants=int(sum(neighbours)),isolated=isolated,
                          smoothing_sigma_arcsec=kernel,psf=psf,psf_matched=psf["psf_matched"],center_offset_arcsec=cube_sep,
                          nexp_center=int(arrays[3][int(math.floor(y+.5)), int(math.floor(x+.5))]),
                          maskbits_center=int(arrays[2][int(math.floor(y+.5)), int(math.floor(x+.5))]),
                          broad_aperture=broad,aperture_radius_arcsec=1.5,
                          eligible=bool(result["status"]=="measured" and psf["psf_matched"] and isolated
                                        and cube_sep is not None and cube_sep<.3
                                        and abs(float(cube.get("aperture_radius_arcsec",0))-1.5)<.05))
            return result
    except Exception as exc:
        return {**fail,**provenance,"reason":str(exc)}


def forced_photometry(rec,legacy,cache_dir,cube=None,*,bands=("r","z")):
    """Return a band-keyed dict. Use the same cache path/object as legacy_brick_cones.

    'eligible' refers only to image quality, isolation, and a model-PSF-matched
    1.5-arcsec aperture at the cube position. Interpret continuum/EW, bandpass,
    foreground correction and line uncertainty in the calling scientific analysis.
    """
    if not all(b in ("g","r","z") for b in bands):
        raise ValueError("unsupported native band")
    return {b:_band(rec,legacy,cache_dir,cube,b) for b in bands}
