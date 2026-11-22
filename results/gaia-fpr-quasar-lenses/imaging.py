"""Pixel-level foreground-light search and SIE + external shear consistency test.

All offsets are east/north arcsec. Position angles below run from east to north.
Two-band residuals and a successful fit are candidate evidence, not a redshift
measurement. Doubles have fewer positional constraints than free model parameters.
"""
import hashlib
import io
import itertools
import json
import math
import re
import time
from email.utils import parsedate_to_datetime
from threading import Lock
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
from PIL import Image, ImageDraw
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import brentq, least_squares, linear_sum_assignment
from scipy.spatial.distance import cdist, pdist


def sky_offsets(coords, origin):
    """Tangent offsets; exact spherical transformation handles RA wrap and poles."""
    # FITS readers may attach an observation time to an ICRS SkyCoord. Astropy
    # rejects spherical_offsets_to across unequal optional metadata even when
    # both sets of RA/Dec values are already fully transformed into ICRS.
    src = coords.icrs
    ref = origin.icrs
    plain_src = SkyCoord(ra=src.ra, dec=src.dec, frame="icrs")
    plain_ref = SkyCoord(ra=ref.ra, dec=ref.dec, frame="icrs")
    lon, lat = plain_ref.spherical_offsets_to(plain_src)
    return np.column_stack((lon.arcsec, lat.arcsec))


def world(origin, xy):
    c = origin.spherical_offsets_by(xy[0] * u.arcsec, xy[1] * u.arcsec)
    return [float(c.ra.deg), float(c.dec.deg)]


def deflection(v, einstein, q, angle):
    """Kormann SIE (area-equivalent Einstein radius). Arrays end in x,y."""
    ca, sa = np.cos(angle), np.sin(angle)
    rot = np.array([[ca, -sa], [sa, ca]])
    x, y = (v @ rot).T
    eps = np.sqrt(1 - q*q)
    psi = np.maximum(np.hypot(q*x, y), 1.e-12)
    norm = einstein * np.sqrt(q) / eps
    return norm * np.column_stack((np.arctan(eps*x/psi), np.arctanh(eps*y/psi))) @ rot.T


def trace(v, p):
    b, q, phi, g1, g2 = p[2:]
    shear = np.array([[g1, g2], [g2, -g1]])
    return v - v @ shear - deflection(v, b, q, phi)


def image_roots(beta, p):
    """For each polar direction a scale-free SIE has radius-independent bending."""
    b, q, phi, g1, g2 = p[2:]
    matrix = np.array([[1-g1, -g2], [-g2, 1+g1]])

    def equation(angles):
        n = np.column_stack((np.cos(angles), np.sin(angles)))
        direction = n @ matrix
        h = beta + deflection(n, b, q, phi)
        cross = direction[:, 0]*h[:, 1] - direction[:, 1]*h[:, 0]
        radius = np.sum(h*direction, axis=1)/np.sum(direction**2, axis=1)
        return cross, radius

    grid = np.linspace(-np.pi, np.pi, 721)
    residual, _ = equation(grid)
    images, magnifications = [], []
    for i in np.flatnonzero(residual[:-1] * residual[1:] <= 0):
        ang = brentq(lambda a: equation([a])[0][0], grid[i], grid[i+1])
        radius = equation([ang])[1][0]
        v = radius * np.array([np.cos(ang), np.sin(ang)])
        if radius <= 1.e-6 or any(np.linalg.norm(v-w) < 1.e-4 for w in images):
            continue
        h = 1.e-5
        jac = np.column_stack([(trace((v+d)[None], p)[0] - trace((v-d)[None], p)[0])/(2*h)
                               for d in np.eye(2)*h])
        det = np.linalg.det(jac)
        images.append(v)
        magnifications.append(float(1. / det) if abs(det) > 1.e-12 else 1.e12)
    return np.asarray(images).reshape(-1, 2), magnifications


def fit_lens(xy, origin, measured_center=None, g_magnitudes=None):
    sep = float(max(pdist(xy)))
    center = np.mean(xy, axis=0) if measured_center is None else measured_center
    center_basis = "two-band extended-light centroid" if measured_center is not None else "image barycentre"
    if measured_center is None and len(xy)==2:
        # An exact midpoint makes an equal double source coincident with the
        # SIS caustic. Choose one member of the SIS double family to centre the
        # weak prior. Equal flux gives little directional information.
        mags = np.asarray(g_magnitudes if g_magnitudes is not None else [np.nan]*2, dtype=float)
        if np.all(np.isfinite(mags)):
            flux = 10**(-.4*(mags-min(mags)))
            contrast = float((flux[1]-flux[0])/sum(flux))
        else:
            contrast = 0.
        direction = 1 if contrast >= 0 else -1
        center = center + (xy[0]-xy[1]) * .5*direction*np.clip(abs(contrast),.10,.75)
        center_basis = "illustrative double family, displaced towards the fainter image using G contrast"
    center_sigma = max(.35*sep, .2) if measured_center is None else .16
    lower = np.r_[center-2.5*center_sigma, .12*sep, .45, -np.pi, -.15, -.15]
    upper = np.r_[center+2.5*center_sigma, 1.3*sep, .9999, np.pi, .15, .15]
    sigma = .035

    def target(p):
        beta_all = trace(xy-p[:2], p)
        # Gaussian penalties regularize the underdetermined ellipticity/shear.
        return np.r_[((beta_all - np.mean(beta_all, axis=0))/sigma).ravel(),
                     (p[:2]-center)/center_sigma, (p[3]-.85)/.20, p[5:]/.065]

    results = []
    # Offset the unmeasured prior centre for starts only to avoid a critical SIS ring.
    for angle in (-1.05, -.18, .65, 1.45):
        start = np.r_[center + .08*min(center_sigma,sep)*np.array([np.cos(angle), np.sin(angle)]),
                      .5*sep, .87, angle, .011, -.009]
        fit = least_squares(target, start, bounds=(lower, upper), max_nfev=180)
        p = fit.x
        source = trace(xy-p[:2], p).mean(axis=0)
        roots, mu = image_roots(source, p)
        roots += p[:2]
        offsets = cdist(xy, roots)
        oi, ri = linear_sum_assignment(offsets)
        distances = offsets[oi, ri]
        largest_error = float(np.max(distances)) if len(distances) else None
        tolerance = min(.15, max(.065, .035*sep))
        missing_bright = []
        if len(ri) == len(xy):
            for i in sorted(set(range(len(roots))) - set(ri)):
                if abs(mu[i]) >= .5*min(abs(mu[j]) for j in ri) and min(offsets[:, i]) > .55:
                    missing_bright.append(i)
        good = bool(len(oi) == len(xy) and largest_error is not None and largest_error <= tolerance
                    and min(mu, default=0) < 0 < max(mu, default=0)
                    and np.hypot(*p[5:]) <= .15 and not missing_bright)
        result = dict(physical_geometry=good, n_images_model=len(roots),
                      observed_to_predicted=list(zip(oi.tolist(), ri.tolist())),
                      maximum_image_residual_arcsec=largest_error, tolerance_arcsec=tolerance,
                      missing_bright_predictions=missing_bright, theta_e_arcsec=float(p[2]),
                      q=float(p[3]), mass_pa_deg_east_to_north=float(np.degrees(p[4]) % 180),
                      shear1=float(p[5]), shear2=float(p[6]), shear=float(np.hypot(*p[5:])),
                      shear_pa_deg_east_to_north=float(.5*np.degrees(np.arctan2(p[6],p[5]))%180),
                      center_radec_deg=world(origin, p[:2]), source_radec_deg=world(origin, p[:2]+source),
                      source_offset_from_mass_center_arcsec=source.tolist(),
                      predicted_radec_deg=[world(origin, r) for r in roots],
                      signed_magnifications=mu, regularized_objective=float(2*fit.cost),
                      light_center_constrained=measured_center is not None, center_prior_basis=center_basis,
                      positional_dof_without_priors=2*len(xy)-9,
                      interpretation="One regularized solution; q and shear are prior-sensitive. "
                                     "Coincident source positions alone do not establish lensing.")
        results.append(result)
    return min(results, key=lambda r: (not r['physical_geometry'], r['regularized_objective']))


PACE_LOCK = Lock()
NEXT_REQUEST = {}
COOLDOWN = {}


def admit(host):
    while True:
        with PACE_LOCK:
            now = time.monotonic()
            denied = COOLDOWN.get(host,0)-now
            if denied>8:
                raise RuntimeError(f"{host}: cooperative backoff active for {denied:.0f}s")
            delay = max(NEXT_REQUEST.get(host,0)-now, denied)
            if delay<=0:
                spacing = 1.6 if host.endswith("legacysurvey.org") else .35
                NEXT_REQUEST[host] = now+spacing
                return
        time.sleep(delay)


def retrieve(url, params, manifest):
    host = urlparse(url).hostname
    for attempt in range(2):
        try:
            admit(host)
            resp = requests.get(url, params=params, timeout=(6, 24),
                                headers={"User-Agent":"Open-Space-Gaia-FPR-study/1"})
            if resp.status_code in (429, 502, 503, 504):
                header = resp.headers.get("Retry-After", "")
                try:
                    wait = float(header)
                except ValueError:
                    try:
                        wait = parsedate_to_datetime(header).timestamp() - time.time()
                    except (ValueError, TypeError, OverflowError):
                        wait = 20 if resp.status_code==429 else 3
                with PACE_LOCK:
                    COOLDOWN[host] = max(COOLDOWN.get(host,0),time.monotonic()+max(1,wait))
                if attempt==0 and wait<=8:
                    continue
            resp.raise_for_status()
            manifest.append(dict(url=resp.url, bytes=len(resp.content),
                                 sha256=hashlib.sha256(resp.content).hexdigest()))
            return resp.content
        except (requests.ConnectionError, requests.Timeout):
            if attempt:
                raise


def cutouts(origin, survey, manifest, errors, size=192):
    coords = dict(ra=round(float(origin.ra.deg), 8), dec=round(float(origin.dec.deg), 8))
    if survey in ("Legacy-DR10", "Legacy-DR9-North"):
        layer = "ls-dr9-north" if survey=="Legacy-DR9-North" else "ls-dr10"
        payload = retrieve("https://www.legacysurvey.org/viewer/fits-cutout",
                           dict(**coords, size=size, pixscale=.262, layer=layer, bands="rz"), manifest)
        with fits.open(io.BytesIO(payload), memmap=False) as hdus:
            hdu = hdus[0]
            arr = np.asarray(hdu.data, dtype=float)
            bands = str(hdu.header.get("BANDS", "rz")).strip().replace(" ", "")
            if arr.ndim != 3 or len(bands) != len(arr) or any(b not in "rz" for b in bands):
                raise ValueError("unexpected Legacy band cube")
            wcs = WCS(hdu.header, naxis=2)
            for band, image in zip(bands, arr):
                yield band, image, wcs, 22.5
    elif survey.startswith("NOIRLab"):
        release = "ls_dr10" if survey.endswith("DR10") else "coadd/ls_dr9"
        extent = size*.262/3600
        params = dict(POS=f"{coords['ra']},{coords['dec']}",
                      SIZE=f"{extent/max(.01,np.cos(origin.dec.radian))},{extent}",VERB=2)
        payload = retrieve("https://datalab.noirlab.edu/sia/"+release, params, manifest)
        rows = Table.read(io.BytesIO(payload), format="votable")
        rows.convert_bytestring_to_unicode()
        seen = set()
        for row in rows:
            bandpass = str(row["obs_bandpass"]).strip().lower()
            if not bandpass:
                continue
            band = bandpass[0]
            kind = str(row["prodtype"]).strip().lower()
            if band not in "rz" or band in seen or kind!="image":
                continue
            if str(row["proctype"]).strip().lower() not in ("stack","coadd"):
                continue
            try:
                url = str(row["access_url"]).replace("http://","https://",1)
                stamp = retrieve(url, {}, manifest)
                with fits.open(io.BytesIO(stamp),memmap=False) as hdus:
                    hdu = next(h for h in hdus if h.data is not None and h.data.ndim==2)
                    wcs = WCS(hdu.header).celestial
                    if max(hdu.data.shape)>2*size or min(hdu.data.shape)<.65*size:
                        raise ValueError(f"unexpected SIA cutout size {hdu.data.shape}")
                    # Coadd BUNIT and MAGZERO are checked rather than assuming
                    # a calibration for arbitrary images returned by the SIA.
                    unit = str(hdu.header.get("BUNIT","")).lower()
                    zp = 22.5 if "nanomag" in unit else None
                    yield band, np.asarray(hdu.data,dtype=float),wcs,zp
                    seen.add(band)
            except Exception as exc:
                errors.append(f"{survey}/{band}: {type(exc).__name__}: {exc}")
        if not seen:
            raise ValueError("no usable r/z Legacy coadded science image in SIA response")
    elif survey=="PS1" and origin.dec.deg >= -31:
        listing = retrieve("https://ps1images.stsci.edu/cgi-bin/ps1filenames.py",
                           dict(**coords, filters="iz", type="stack"), manifest)
        table = Table.read(listing.decode("ascii").splitlines(), format="ascii")
        for band in "iz":
            for row in table[table["filter"] == band][:1]:
                try:
                    payload = retrieve("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi",
                                       dict(**coords, size=size, format="fits", red=str(row["filename"])), manifest)
                    with fits.open(io.BytesIO(payload), memmap=False) as hdus:
                        hdu = hdus[0]
                        # MAST PS1 "Stack Images" and "Image Cutout Service":
                        # fitscut serves linear fluxes, unlike full asinh skycells.
                        t = hdu.header.get("EXPTIME", 0)
                        zeropoint = 25+2.5*np.log10(t) if t > 0 else None
                        yield band, np.asarray(hdu.data, dtype=float), WCS(hdu.header).celestial, zeropoint
                except Exception as exc:
                    errors.append(f"{survey}/{band}: {type(exc).__name__}: {exc}")


def mad(data):
    return 1.4826 * np.nanmedian(np.abs(data - np.nanmedian(data)))


def moffat(r2, fwhm):
    # beta=3; continuous unit-flux profile; sampling is adequate at >2 pixels/FWHM.
    alpha2 = fwhm**2/(4*(2**(1/3)-1))
    return 2/(np.pi*alpha2) * (1+r2/alpha2)**(-3)


def field_psf(image, stars_xy, sky, noise, scale):
    # Use compact isolated off-target peaks to bracket the field PSF; galaxies
    # broaden the distribution, hence the lower quartile instead of the median.
    signal = np.nan_to_num((image-sky)/noise, nan=-1.e6)
    peaks = (signal == maximum_filter(signal, size=17)) & (signal > 18)
    yy, xx = np.nonzero(peaks)
    widths = []
    y, x = np.mgrid[-7:8, -7:8]
    for i in np.argsort(signal[yy,xx])[::-1][:30]:
        cy,cx = yy[i],xx[i]
        if min(cy,cx,image.shape[0]-1-cy,image.shape[1]-1-cx) < 8 or min(
                np.linalg.norm(stars_xy-[cx,cy], axis=1)) < 9/scale:
            continue
        patch = signal[cy-7:cy+8, cx-7:cx+8]
        if not np.all(np.isfinite(image[cy-7:cy+8, cx-7:cx+8])):
            continue
        amp = signal[cy,cx]
        def fun(p):
            model = p[0] * (1+4*(2**(1/3)-1)*((x-p[1])**2+(y-p[2])**2)/p[3]**2)**(-3)+p[4]
            return ((model-patch)/(1+.035*np.maximum(patch,0))).ravel()
        fit = least_squares(fun, [amp,0,0,1.2/scale,0],
                            bounds=([.3*amp,-1,-1,.50/scale,-5],[2*amp,1,1,2.4/scale,5]), max_nfev=30)
        if np.mean(fun(fit.x)**2) < 3.5 and .53 < fit.x[3]*scale < 2.3:
            widths.append(fit.x[3])
    return (float(np.percentile(widths, 25)) if widths else None), len(widths)


def band_fit(image, wcs, coords, origin, zeropoint):
    scale = float(np.mean(proj_plane_pixel_scales(wcs))*3600)
    if not .20 < scale < .31 or image.ndim != 2:
        raise ValueError("unexpected pixel scale or FITS dimensions")
    pixels = np.column_stack(wcs.world_to_pixel(coords))
    cx, cy = np.mean(pixels, axis=0)
    sep = float(np.max(pdist(pixels)))
    h = int(np.ceil(sep/2+3.3/scale))
    xc, yc = int(round(cx)), int(round(cy))
    if min(xc, yc, image.shape[1]-1-xc, image.shape[0]-1-yc) < h+2:
        raise ValueError("group exceeds stamp")
    yy, xx = np.indices(image.shape)
    sky_pix = image[(np.hypot(xx-cx, yy-cy) > h+3) & np.isfinite(image)]
    if len(sky_pix) < 300 or not np.isfinite(sky := np.median(sky_pix)) or not (noise := mad(sky_pix)) > 0:
        raise ValueError("missing or uniform sky")
    field_width, num_stars = field_psf(image, pixels, sky, noise, scale)
    patch = (image[yc-h:yc+h+1, xc-h:xc+h+1] - sky) / noise
    yy, xx = np.indices(patch.shape)
    xx, yy = xx-h, yy-h
    positions = pixels - [xc, yc]
    valid = np.isfinite(patch) & (np.hypot(xx,yy) <= h)
    if valid.sum() < .90*np.pi*h*h or any(np.mean(np.isfinite(patch[
            (xx-x)**2+(yy-y)**2 < (1/scale)**2])) < .95 for x,y in positions):
        raise ValueError("masked or invalid image core")
    # Hybrid empirical sky + fractional core error permits PSF systematics.
    uncertainty = np.hypot(1., .035*np.maximum(np.nan_to_num(patch), 0))
    n = len(positions)
    fmin, fmax = (.64/scale, 2.15/scale) if field_width is None else (.85*field_width, 1.15*field_width)
    finit = 1.2/scale if field_width is None else field_width
    amp_guess = [max(12, patch[int(round(y))+h,int(round(x))+h] /
                     moffat(0, finit)) for x,y in positions]
    broad = max(float(np.max(patch[valid])), 1)*np.pi*(2*fmax)**2
    highs = [max(100,4*broad)]*n
    low = [0]*n + [-.6/scale, -.6/scale, fmin, -5, -.4, -.4]
    high = highs + [.6/scale,.6/scale,fmax,5,.4,.4]
    initial = np.r_[np.minimum(amp_guess, np.array(highs)*.8), 0,0,finit,0,0,0]
    gx0, gy0 = positions.mean(axis=0)
    window = max(.35/scale, .38*sep)

    def solve(kind, p0, lo, hi):
        def model(p):
            dx,dy,f,b,bx,by = p[n:n+6]
            arr = b + bx*xx + by*yy
            for flux,(x,y) in zip(p[:n],positions):
                arr += flux*moffat((xx-x-dx)**2+(yy-y-dy)**2, f)
            if kind is not None:
                if isinstance(kind,int):
                    flux = p[n+6]
                    gx,gy = positions[kind] + [dx,dy]
                else:
                    flux,gx,gy = p[n+6:n+9]
                g = moffat((xx-gx)**2+(yy-gy)**2, f)
                if kind=="extended" or isinstance(kind,int):
                    g = gaussian_filter(g, p[-1], mode="constant", truncate=3.5)
                arr += flux*g
            return arr

        def residual(p):
            return ((model(p)-patch)/uncertainty)[valid]
        fit = least_squares(residual, p0, bounds=(lo,hi), max_nfev=80,
                            ftol=1.e-5, xtol=1.e-5, x_scale="jac")
        return fit, model(fit.x)

    base, base_img = solve(None, initial, low, high)
    init = np.r_[base.x, max(5, min(broad*.07,np.sum(base.x[:n])*.18)),gx0,gy0]
    lo = low+[0, gx0-window, gy0-window]
    hi = high+[4*broad, gx0+window, gy0+window]
    point, _ = solve("point", init, lo, hi)
    extended, ext_img = solve("extended", np.r_[point.x, .55/scale],
                              lo+[.12/scale], hi+[2.0/scale])
    p = extended.x
    gx, gy = p[n+7:n+9]
    # Remove the fitted common astrometric shift before comparing with Gaia.
    gal_px = np.array([gx+xc-p[n],gy+yc-p[n+1]])
    gal = wcs.pixel_to_world(*gal_px).icrs
    near = float(min(np.linalg.norm(pixels-gal_px, axis=1))*scale)
    center_offset = float(np.linalg.norm(gal_px-np.mean(pixels, axis=0))*scale)
    # Estimate noise correlation via residuals of same weight in blank sky.
    # Conservative floor of 2.5 on variance inflation from resampling.
    blank = (np.hypot(*np.meshgrid(np.arange(image.shape[1])-cx, np.arange(image.shape[0])-cy)) > h+3)
    signal = (image-sky)/noise
    blank &= np.isfinite(signal) & (np.abs(signal)<3)
    variance_inflation = 1.
    for dy,dx in [(dy,dx) for dy in range(0,4) for dx in range(-3,4) if dy>0 or dx>0]:
        rolled = np.roll(signal,(dy,dx),axis=(0,1))
        ok = blank & np.roll(blank,(dy,dx),axis=(0,1))
        ok[:dy,:] = False
        if dx>0:
            ok[:,:dx] = False
        elif dx<0:
            ok[:,dx:] = False
        if ok.sum()>100:
            variance_inflation += 2*max(0, np.mean(signal[ok]*rolled[ok])/np.mean(signal[ok]**2))
    variance_inflation = max(2.5, float(variance_inflation))
    penalty = np.log(max(valid.sum()/variance_inflation, 2))
    improvement = 2*(base.cost-extended.cost)/variance_inflation - 4*penalty
    third_vs_extended = 2*(point.cost-extended.cost)/variance_inflation - penalty
    covariance = np.linalg.pinv(extended.jac.T @ extended.jac, rcond=1.e-12)
    errors = np.sqrt(np.maximum(np.diag(covariance),0)*variance_inflation*
                     max(1, 2*extended.cost/(valid.sum()-len(p))))
    snr = np.divide(p[:n],errors[:n],out=np.zeros(n),where=errors[:n]>0)
    gal_snr = float(p[n+6]/errors[n+6]) if errors[n+6]>0 else 0
    size_arcsec = float(p[-1]*scale)
    legible = bool(num_stars > 0 and min(snr) > 4 and
                   (np.max(np.abs((patch-ext_img)[valid]/uncertainty[valid])) < 13)
                   and (2*extended.cost/valid.sum() < 6))
    evidence = bool(legible and gal_snr >= 5 and improvement >= 14 and third_vs_extended >= 8
                    and near >= max(.35,.35*p[n+2]*scale) and size_arcsec >= .27
                    and size_arcsec < 1.95 and center_offset <= .45*sep*scale)
    host_scores = []
    if evidence:
        # A physical quasar pair can have extended host light at either image.
        # Compare with the same PSF and convolved profile, locking only the
        # extended component's centroid to one shifted Gaia image in each fit.
        for i in range(n):
            host, _ = solve(i,np.r_[p[:n+7],p[-1]],low+[0,.12/scale],high+[4*broad,2./scale])
            score = 2*(host.cost-extended.cost)/variance_inflation - 2*penalty
            host_scores.append(float(score))
    mags = [(float(zeropoint-2.5*np.log10(flux*noise)) if zeropoint is not None and flux>0 else None)
            for flux in p[:n]]
    diag = dict(psf_search_width_arcsec=image.shape[1]*scale, field_psf_objects=num_stars, empirical_fwhm_arcsec=None if field_width is None
                else float(field_width*scale), fitted_fwhm_arcsec=float(p[n+2]*scale),
                point_magnitudes_AB=mags, point_snr=snr.tolist(),
                point_flux_in_native_units=(p[:n]*noise).tolist(),
                pixel_scale_arcsec=scale, background_sky_per_pixel=float(sky),
                sky_rms_per_pixel=float(noise), variance_inflation=variance_inflation,
                bic_points_minus_extended=float(improvement), bic_third_point_minus_extended=float(third_vs_extended),
                galaxy_snr=gal_snr, galaxy_radec_deg=[float(gal.ra.deg), float(gal.dec.deg)],
                galaxy_offset_arcsec=sky_offsets(SkyCoord([gal]), origin)[0].tolist(),
                galaxy_intrinsic_gaussian_sigma_arcsec=size_arcsec, nearest_image_arcsec=near,
                bic_single_host_minus_free_galaxy=host_scores,
                resolved_from_every_single_host=bool(host_scores and min(host_scores)>=8),
                usable_morphology=legible, extended_light=evidence)
    return diag, [patch, patch-base_img, patch-ext_img]


def save_panel(path, panels):
    # Top-to-bottom filters; left-to-right image, residual from points, residual with galaxy.
    thumb = 128
    canvas = Image.new("L", (3*thumb, (thumb+18)*len(panels)), 22)
    draw = ImageDraw.Draw(canvas)
    for i,(label, arrays) in enumerate(panels):
        y = (thumb+18)*i
        draw.text((2,y+2),label+"     data / points residual / galaxy residual",fill=240)
        for j,a in enumerate(arrays):
            arr = np.nan_to_num(a, nan=0.)
            lim = max(5.,np.percentile(np.abs(arr), 99))
            gray = np.clip(127+120*np.arcsinh(arr/2)/np.arcsinh(lim/2),0,255).astype(np.uint8)
            canvas.paste(Image.fromarray(gray[::-1]).resize((thumb,thumb)), (j*thumb,y+18))
    canvas.save(path)


def assess_one(system, output_dir):
    ident = str(system["id"])
    coords = SkyCoord([float(c["ra"]) for c in system["components"]],
                      [float(c["dec"]) for c in system["components"]], unit="deg")
    if len(coords)<2:
        raise ValueError("at least two image positions required")
    shift = sky_offsets(coords, coords[0]).mean(axis=0)
    origin = SkyCoord(*world(coords[0], shift), unit="deg")
    xy = sky_offsets(coords, origin)
    pairs = [dict(i=i,j=j,arcsec=float(coords[i].separation(coords[j]).arcsec))
             for i,j in itertools.combinations(range(len(coords)),2)]
    records, errors, stamps, manifest = {}, [], {}, []
    north = origin.dec.deg >= 32.375 and origin.galactic.b.deg > 0
    viewer = "Legacy-DR9-North" if north else "Legacy-DR10"
    alternatives = ("NOIRLab-DR9",) if north else ("NOIRLab-DR10","NOIRLab-DR9")
    for survey in (viewer,"PS1",*alternatives):
        try:
            for band, img, wcs, zp in cutouts(origin, survey, manifest, errors):
                key = survey+"/"+band
                try:
                    records[key], stamps[key] = band_fit(img,wcs,coords,origin,zp)
                except Exception as exc:
                    errors.append(f"{key}: {type(exc).__name__}: {exc}")
        except Exception as exc:
            errors.append(f"{survey}: {type(exc).__name__}: {exc}")
        need_psf = {k for k,v in records.items() if k.startswith(survey+"/") and v["field_psf_objects"]==0}
        if need_psf:
            try:
                # A 120-arcsec field samples 6.25 times the calibration area.
                for band, img, wcs, zp in cutouts(origin, survey, manifest, errors, size=480):
                    key = survey+"/"+band
                    if key in need_psf:
                        try:
                            wider, panel = band_fit(img, wcs, coords, origin, zp)
                            if wider["field_psf_objects"]:
                                records[key], stamps[key] = wider, panel
                        except Exception as exc:
                            errors.append(f"{key} wide: {type(exc).__name__}: {exc}")
            except Exception as exc:
                errors.append(f"{survey} wide: {type(exc).__name__}: {exc}")
        if len({k.split("/")[-1] for k,x in records.items() if x["usable_morphology"]}) >= 2:
            break

    good = sorted([k for k,x in records.items() if x["extended_light"]])
    usable = [k for k,x in records.items() if x["usable_morphology"]]
    photometric_tests = []
    for a,b in itertools.combinations(usable,2):
        if a.split("/")[-1]==b.split("/")[-1]:
            continue
        f1 = np.asarray(records[a]["point_flux_in_native_units"])
        f2 = np.asarray(records[b]["point_flux_in_native_units"])
        spread = float(np.ptp(-2.5*np.log10(f1/f2)))
        photometric_tests.append(dict(bands=[a,b], colour_spread_mag=spread,
                                      within_055mag=spread<=.55))
    selected, chromaticity = None, None
    for a,b in itertools.combinations(good, 2):
        if a.split("/")[-1]==b.split("/")[-1]:
            continue
        first, second = records[a],records[b]
        d = np.linalg.norm(np.subtract(first["galaxy_offset_arcsec"],second["galaxy_offset_arcsec"]))
        if d>.35:
            continue
        # Colour spread among the Gaia images, independent of absolute zeropoints.
        color_difference = -2.5*np.log10(np.array(first["point_flux_in_native_units"]) /
                                         np.array(second["point_flux_in_native_units"]))
        color_spread = float(np.ptp(color_difference))
        if color_spread<=.55 and (selected is None or color_spread<chromaticity):
            selected, chromaticity = (a,b), color_spread
    evidence_center = None
    if selected:
        evidence_center = np.mean([records[k]["galaxy_offset_arcsec"] for k in selected],axis=0)
        status = "extended_light"
    elif not records:
        status = "unavailable"
    elif len({k.split("/")[-1] for k in usable}) < 2 or good:
        status = "inconclusive"
    else:
        status = "no_robust_extended_light"
    try:
        fit = fit_lens(xy, origin, evidence_center, [c.get("g") for c in system["components"]])
    except Exception as exc:
        errors.append(f"lens model: {type(exc).__name__}: {exc}")
        fit = {"physical_geometry": False, "evaluation_failed": True}
    result = dict(id=ident, imaging_status=status, imaging_model_pass=bool(selected and fit["physical_geometry"]),
                  pair_separations=pairs, bands=records, matched_bands=selected,
                  image_colour_spread_mag=chromaticity, lens_model=fit, imaging_inputs=manifest,
                  two_band_photometry=photometric_tests,
                  consistent_two_band_photometry=any(t["within_055mag"] for t in photometric_tests),
                  free_centroid_preferred_to_each_single_host_in_both_bands=bool(selected and all(
                      records[k]["resolved_from_every_single_host"] for k in selected)),
                  failures=errors)
    if selected or good:
        filename = re.sub(r"[^A-Za-z0-9_.-]", "_", ident)[:60]+"-"+hashlib.sha256(ident.encode()).hexdigest()[:8]+".png"
        save_panel(output_dir/filename, [(key, stamps[key]) for key in records])
        result["diagnostic_png"] = filename
    return result


def assess(systems, output_dir, *, workers=6):
    """Run identical fitting on the supplied target/control/positive-control stream."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    def scrub(x):
        if isinstance(x, dict):
            return {k:scrub(v) for k,v in x.items()}
        if isinstance(x, (list, tuple)):
            return [scrub(v) for v in x]
        if isinstance(x, (float, np.floating)):
            return float(x) if np.isfinite(x) else None
        if isinstance(x, np.generic):
            return x.item()
        return x

    def run(system):
        try:
            return scrub(assess_one(system, output_dir))
        except Exception as exc:
            return dict(id=str(system["id"]), imaging_status="inconclusive", imaging_model_pass=False,
                        lens_model={"physical_geometry":False}, failures=[f"{type(exc).__name__}: {exc}"])
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as pool, (output_dir/"imaging.jsonl").open("w") as f:
        for item in pool.map(run, systems):
            f.write(json.dumps(item, allow_nan=False)+"\n")
            f.flush()
            rows.append(item)
    return rows
