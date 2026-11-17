"""RDI, point-source search and injection recovery on public Ross 154 F2100W mosaics.

python ross_rdi.py FITS_CACHE OUTPUT_DIR. Verified bytes from inputs.json are
downloaded as needed. Requires numpy, scipy, astropy. Coordinates are east/north.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi
from scipy.optimize import least_squares
from scipy.signal import correlate
from scipy.stats import t

HERE = Path(__file__).resolve().parent
HALF = 106
FWHM = 6.3  # approx. 0.70 arcsec for F2100W
RADII = [3, 4.5, 6, 7.5, 9]  # edge in arcsec
PA = np.arange(48) * (2 * np.pi / 48)
yy, xx = np.mgrid[-HALF:HALF+1, -HALF:HALF+1]
radius = np.hypot(yy, xx)


def write_csv(path, records):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def fits_for(item, cache):
    p = cache / item["filename"]
    if not p.exists():
        with urlopen(item["url"], timeout=90) as resp:
            p.write_bytes(resp.read())
    if hashlib.sha256(p.read_bytes()).hexdigest() != item["sha256"]:
        raise ValueError(f"Input changed: {p}")
    return p


def frame(item, cache):
    with fits.open(fits_for(item, cache)) as f:
        img = f['SCI'].data.astype(float)
        wht = f['WHT'].data.astype(float)
        h = f['SCI'].header.copy()
        mjd = h['MJD-AVG']
        skyval = np.nanmedian(img[wht > 0.6 * np.nanmax(wht)])
        img -= skyval
        cy, cx = np.unravel_index(np.nanargmax(img), img.shape)
        dy, dx = np.mgrid[-5:6,-5:6]
        cut = img[cy-5:cy+6,cx-5:cx+6]
        def resid(par):
            bg, amp, x, y, sigma = par
            return (bg + amp * np.exp(-((dx-x)**2 + (dy-y)**2)/(2*sigma**2)) - cut).ravel()
        fit = least_squares(resid, [0, cut.max(), 0, 0, 2.7]).x
        cx, cy = cx + fit[2], cy + fit[3]
        sample = (yy + cy, xx + cx)
        valid = ndi.map_coordinates(np.isfinite(img).astype(float) * (wht > 0.30 * np.nanmax(wht)),
                                    sample, order=1, mode="constant") > .99
        arr = ndi.map_coordinates(np.nan_to_num(img), sample, order=1, mode="constant")
        w = WCS(h, fix=False)
        mat = w.pixel_scale_matrix * 3600
        return dict(image=arr, valid=valid, cx=cx, cy=cy, mjd=mjd, mat=mat,
                    scale=abs(np.linalg.det(mat))**.5, jwst_xyz=np.array([h['JWST_'+k] for k in 'XYZ']),
                    center_adu=float(img[int(round(cy)),int(round(cx))]),
                    weight_core_ratio=wht[int(round(cy)),int(round(cx))]/np.nanmax(wht))


def reduce(sci, ref, valid, psf):
    # Robust host normalization uses only the core and first Airy ring, outside all injection sites.
    zone = valid & (radius < 19)
    def diff(par):
        return (sci[zone] - par[0] * ref[zone] - par[1])
    p = least_squares(diff, [1.6, 0.], loss="soft_l1", f_scale=8).x
    raw = sci - p[0]*ref - p[1]
    raw[~valid] = 0
    # Zero-sum filter simultaneously fits a constant on a compact core aperture.
    return correlate(raw, psf, mode="same", method="fft") / np.sum(psf**2)


def get_psf(star):
    r = radius[HALF-10:HALF+11, HALF-10:HALF+11]
    profile = star[HALF-10:HALF+11, HALF-10:HALF+11].copy()
    profile -= np.median(profile[(r>=9) & (r<=10)])
    profile[r>10]=0
    # Empirical injection core shape relative to the measured host.
    kernel = profile.copy()
    kernel[r<=10] -= np.mean(kernel[r<=10])
    return profile, kernel


def unit_sig(scale, valid, sigmap):
    sigs = []
    for lo, hi in zip(RADII, RADII[1:]):
        reg = valid & (radius * scale >= lo) & (radius * scale < hi)
        data = sigmap[reg]
        med = np.median(data)
        sigma = 1.4826 * np.median(abs(data-med))
        # Independent elements approx. one FWHM on a circle. One-sided Gaussian 5-sigma tail.
        fraction = reg.sum() / ((radius * scale >= lo) & (radius * scale < hi)).sum()
        n = max(3, int(fraction * 2 * np.pi * lo / (FWHM*scale)))
        th = t.ppf(1-2.86651572e-7, n-1)*np.sqrt(1+1/n)
        sigs.append((reg, med, sigma, th))
    return sigs


def peaks(sigmap, sigs):
    found = []
    local = sigmap == ndi.maximum_filter(sigmap, size=7)
    for reg, med, sigma, th in sigs:
        for iy, ix in np.argwhere(reg & local & (sigmap >= med + sigma*th)):
            # Non-maximum suppression prevents duplicates for the same resolution element.
            found.append((float((sigmap[iy,ix]-med)/sigma), ix, iy, sigmap[iy,ix]))
    keep = []
    for rec in sorted(found, reverse=True):
        if all(np.hypot(rec[1]-x, rec[2]-y)>FWHM for _,x,y,_ in keep):
            keep.append(rec)
    return keep


def centroid_peak(sigmap, ix, iy):
    # Quadratic interpolation along each axis near the maximum.
    def off(a,b,c):
        den = a-2*b+c
        return np.clip((a-c)/(2*den),-.5,.5) if den<0 else 0.
    return ix+off(*sigmap[iy,ix-1:ix+2]), iy+off(*sigmap[iy-1:iy+2,ix])


def main(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    items = json.loads((HERE/'inputs.json').read_text())
    frames = {row["obs_id"]: frame(row,cache) for row in items}
    refrows = [frames[r['obs_id']] for r in items if r['target_name']=="Wolf-359"]
    ref = np.mean([r['image'] for r in refrows], axis=0)
    refvalid = np.logical_and.reduce([r['valid'] for r in refrows])
    curve, sources, meta, trials, masks, qa = [], [], [], [], [], []
    sky_n, sky_e = np.mgrid[-9:9.001:.05,-9:9.001:.05]
    sky_pts = np.array([sky_e.ravel(), sky_n.ravel()])
    with (HERE/'catalogue/gaia_dr3.csv').open() as f:
        host = next(r for r in csv.DictReader(f) if r['source_id']=='4075141768785646848')
    for row in sorted(items, key=lambda r:r['t_min']):
        if row['target_name'] != "Ross-154":
            continue
        name = row['obs_id'].split('-')[1].split('_')[0]
        f=frames[row['obs_id']]
        img, scale = f['image'], f['scale']
        valid = refvalid & f['valid']
        core, kernel = get_psf(img)
        data = reduce(img,ref,valid,kernel)
        # Exclude any footprint location without fully exposed support for the filter.
        support = ndi.minimum_filter(valid.astype(int), size=21).astype(bool)
        ix,iy = np.rint(np.linalg.solve(f['mat'], sky_pts)+HALF).astype(int)
        masks.append(support[iy,ix])
        sigs = unit_sig(scale,support,data)
        peaklist = peaks(data,sigs)
        for snr,ix,iy,contrast in peaklist:
            x,y=centroid_peak(data,ix,iy)
            east,north=f['mat'] @ np.array([x-HALF,y-HALF])
            sources.append(dict(obs=name,id=len(sources)+1,east_as=round(east,4),north_as=round(north,4),
                                detection_snr=round(snr,2),core_contrast=round(contrast,7)))
        # The RDI scaling domain does not overlap injections at r >=3 arcsec (minimum
        # injection radius 3.75 arcsec, footprint radius 10 pixels). Reduction after
        # RDI is linear. An injected source changes the filtered image by the translated
        # core/kernel cross-correlation exactly. Work at integer output-pixel locations.
        reference_core,_ = get_psf(refrows[1]['image'])
        # Same finite aperture convention as primary host photometry. Inject a different
        # star's empirical PSF; the detection kernel uses the target's measured shape.
        reference_core *= core.sum()/reference_core.sum()
        for j,(reg,med,sigma,th) in enumerate(sigs):
            lo, hi = RADII[j:j+2]
            rad = (lo+hi)/2/scale
            # Uniform position angles fixed in sky; no masks around genuine dark/bright sources.
            sites=[]
            for a,pa in enumerate(PA + 0.043):
                delta=np.linalg.solve(f['mat'],rad*scale*np.array([np.sin(pa),np.cos(pa)]))
                x,y=np.rint(np.array([HALF,HALF])+delta).astype(int)
                if reg[y,x]:
                    dx,dy = np.array([HALF,HALF])+delta - [x,y]
                    stamp = ndi.shift(np.pad(reference_core,2),(dy,dx),order=3,mode='constant')
                    response = correlate(stamp,kernel,mode='full')/np.sum(kernel**2)
                    h=response.shape[0]//2
                    sites.append((x,y,a,stamp,response[h-7:h+8,h-7:h+8]))
            if sites:
                x,y,_,stamp,patch=sites[0]
                injected = img.copy()
                injected[y-12:y+13,x-12:x+13] += th*sigma*stamp
                full=reduce(injected,ref,valid,kernel)
                maxerr=np.max(abs(full[y-7:y+8,x-7:x+8]-data[y-7:y+8,x-7:x+8]-th*sigma*patch))
                assert maxerr < 1e-8, (name, lo, maxerr)
                qa.append(dict(obs=name, radius_arcsec=(lo+hi)/2, max_full_vs_linear_error=float(maxerr)))
            for factor in (.75, 1., 1.25, 1.5, 2.):
                contrast = factor*sigma*th
                rec=0
                for x,y,a,stamp,response in sites:
                    # Matching requires a local maximum within half a FWHM above threshold
                    # that was not already a threshold detection before injection.
                    orig=data[y-7:y+8,x-7:x+8]
                    inj=orig+contrast*response
                    oy,ox=np.mgrid[-7:8,-7:8]
                    area=np.hypot(ox,oy)<=FWHM/2
                    def detected(patch):
                        return bool(np.any(area & (patch==ndi.maximum_filter(patch,7)) &
                                           (patch>med+sigma*th)))
                    recovered=detected(inj) and not detected(orig)
                    rec += recovered
                    trials.append(dict(obs=name,radius_arcsec=(lo+hi)/2,pa_deg=round(np.rad2deg(PA[a]+.043),4),
                                       contrast=round(contrast,8),recovered=int(recovered)))
                curve.append(dict(obs=name,radius_arcsec=round(rad*scale,3),bin_inner_arcsec=lo,
                                  bin_outer_arcsec=hi,coverage_fraction=round(reg.sum()/((radius*scale>=lo)&(radius*scale<hi)).sum(),4),
                                  cutoff_sigma=round(th,3),threshold_contrast=round(th*sigma,8),
                                  injected_contrast=round(contrast,8),recovered=rec,injected=len(sites)))
        # Stationary-background reflex vector at epoch relative to the first science product.
        if not meta:
            first = f
        ra,dec=np.deg2rad([float(host['ra']),float(host['dec'])])
        axes=np.array([[-np.sin(ra),np.cos(ra),0.],
                       [-np.sin(dec)*np.cos(ra),-np.sin(dec)*np.sin(ra),np.cos(dec)]])
        pm=np.array([float(host['pmra']),float(host['pmdec'])])
        bg=-pm*(f['mjd']-first['mjd'])/365.25 + float(host['parallax'])*axes@((f['jwst_xyz']-first['jwst_xyz'])/149597870.7)
        meta.append(dict(obs=name,mjd=f['mjd'],center_x=f['cx'],center_y=f['cy'],pixel_scale=scale,
                         core_weight_fraction=f['weight_core_ratio'],static_bg_east_mas=bg[0],static_bg_north_mas=bg[1],
                         n_peaks=len(peaklist)))
    coverage = []
    sky_r=np.hypot(sky_pts[0],sky_pts[1])
    first,aug1,aug2=masks
    for lo,hi in zip(RADII,RADII[1:]):
        use=(sky_r>=lo)&(sky_r<hi)
        coverage.append(dict(inner_as=lo,outer_as=hi,any_visit=float(np.mean((first|aug1|aug2)[use])),
                             both_epochs=float(np.mean((first & (aug1 | aug2))[use]))))
    write_csv(out/'coverage.csv',coverage)
    write_csv(out/'linearity_check.csv',qa)
    write_csv(out/'injection_trials.csv',trials)
    write_csv(out/'contrast_recovery.csv',curve)
    write_csv(out/'visits.csv',meta)
    if sources:
        write_csv(out/'peaks.csv',sources)
    else:
        (out/'peaks.csv').write_text('obs,id,east_as,north_as,detection_snr,core_contrast\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument("cache",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    main(a.cache,a.out)
