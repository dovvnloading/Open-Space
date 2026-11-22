"""Empirical line and bridge tests on masked, continuum-subtracted PDR1 IFU cubes."""
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

PDR = "https://web.corral.tacc.utexas.edu/hetdex/HETDEX/pdr/pdr1/datacubes/"
MPE = "https://data.sciserver.mpe.mpg.de/public/hetdex/pdr/pdr1/datacubes/"


def mad(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a)))) if len(a) else np.nan


def empirical(value, controls):
    controls = np.array(controls, float)
    controls = controls[np.isfinite(controls)]
    s = mad(controls)
    if len(controls) < 12 or not s > 0 or not np.isfinite(value):
        return dict(n=len(controls), sn=None, value=None, scatter=None, tail_fraction=None)
    return dict(value=float(value - np.median(controls)), scatter=s, n=len(controls),
                sn=float((value-np.median(controls))/s),
                tail_fraction=float((1+np.sum(controls >= value))/(1+len(controls))))


def cube_test(rec, fetch):
    shot, slot = str(int(rec["shotid"])), str(int(rec["ifuslot"])).zfill(3)
    name = f"dex_cube_{shot}_{slot}.fits"
    suffix = f"{shot}/{name}"
    file = fetch("cube", name, [PDR + suffix, MPE + suffix])
    # FITS tile-compressed arrays are read only in the required windows.
    with fits.open(file, memmap=False) as h:
        return _cube_measure(h, rec, name)


def _cube_measure(h, rec, name, halo_trials=False, resident=False):
    hd = h["DATA"]
    wave = np.arange(hd.shape[0]) * 2 + 3470.
    pi = int(np.abs(wave-rec["wave"]).argmin())
    wcs = WCS(hd.header).celestial
    ra = rec.get("ra_det") if rec.get("ra_det") is not None else rec["ra"]
    dec = rec.get("dec_det") if rec.get("dec_det") is not None else rec["dec"]
    x, y = map(float, wcs.world_to_pixel_values(ra, dec))
    scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.)
    seeing = rec.get("fwhm") or hd.header.get("FWHM") or hd.header.get("SEEING")
    if seeing is None or not 0.5 < float(seeing) < 5:
        return {"status":"unmeasured seeing", "cube":name}
    seeing = float(seeing)
    yy, xx = np.mgrid[:hd.shape[1], :hd.shape[2]]
    sep1 = np.hypot(xx-x, yy-y)*scale
    if not np.any(sep1 < 0.75):
        return {"status":"off cube", "cube":name}
    ap = sep1 < 1.5
    outer = (sep1 >= 9) & (sep1 <= 20)
    centers = [(x, y)]
    if rec.get("secondary_ra") is not None:
        centers.append(tuple(map(float, wcs.world_to_pixel_values(
            rec["secondary_ra"], rec["secondary_dec"]))))
    psfs = []
    beta = 2.5
    alpha = seeing / (2*np.sqrt(2**(1/beta)-1))
    for cx, cy in centers:
        r = np.hypot(xx-cx, yy-cy)*scale
        psfs.append((1+(r/alpha)**2)**(-beta))
    core = psfs[0] > 0.05
    if len(centers)==2:
        cx,cy = centers[1]
        sep2 = np.hypot(xx-cx,yy-cy)*scale
        direction = np.array([cx-x,cy-y])
        length = np.hypot(*direction)
        along = ((xx-x)*direction[0]+(yy-y)*direction[1])/length
        across = ((xx-x)*direction[1]-(yy-y)*direction[0])/length
        bridge = ((np.abs(across)*scale < 0.75*seeing) & (along>0) & (along<length)
                  & (sep1>1.3*seeing) & (sep2>1.3*seeing))
        core |= psfs[1] > 0.05
        outer &= sep2>9
    half = 3  # seven native bins, or 14 A
    # Disjoint line windows on both sides, skipping 40 A around the target.
    controls = [pi+j for j in range(-180,181,9) if abs(j)>=23 and 14<=pi+j<len(wave)-14]
    if not 14<=pi<len(wave)-14:
        return {"status":"spectral edge", "cube":name}
    pmin = min(controls+[pi])-13
    pmax = max(controls+[pi])+14
    # Cropping in wavelength also bounds the in-memory decompression.
    def window(hdu):
        view = hdu.data if resident else hdu.section
        return view[pmin:pmax,:,:]

    a = np.array(window(hd), dtype=np.float32, copy=True)
    flags = h["MASK"] if "MASK" in h else h["BITMASK"]
    mask = window(flags)
    err = window(h["ERROR"])
    a[(mask!=0) | ~(err>0) | ~np.isfinite(err) | ~np.isfinite(a)] = np.nan

    def narrowband(k):
        k -= pmin
        line = a[k-half:k+half+1]
        side = np.concatenate((a[k-13:k-6],a[k+7:k+14]))
        good = (np.sum(np.isfinite(line),axis=0)==7) & (np.sum(np.isfinite(side),axis=0)>=10)
        m = np.sum(np.nan_to_num(line),axis=0)-7*np.nanmedian(side,axis=0)
        return m, good & np.isfinite(m)

    region = bridge if len(centers)==2 else (sep1>1.5*seeing) & (sep1<3.0*seeing)
    template = np.stack(psfs,axis=-1)

    def measure(m,good):
        if (good & outer).sum()<100:
            return np.nan, np.nan, np.nan
        sky = float(np.median(m[good & outer]))
        if (good & ap).sum() < 0.90*ap.sum():
            return np.nan, np.nan, sky
        m = m - sky
        central = np.mean(m[ap & good])*ap.sum()
        if region.sum()<8 or (good & region).sum()<0.9*region.sum():
            return central, np.nan, sky
        use = core & good
        if use.sum()<12:
            return central, np.nan, sky
        weights = np.linalg.lstsq(template[use],m[use],rcond=None)[0]
        residual = m - template@weights
        return central, float(np.mean(residual[good & region])*region.sum()), sky

    m, good = narrowband(pi)
    target = measure(m,good)
    offset = np.array([measure(*narrowband(k)) for k in controls])
    centre = empirical(target[0], offset[:,0])
    bridge_stats = empirical(target[1], offset[:,1])
    background = empirical(target[2],offset[:,2])
    result = dict(status="measured", cube=name, central_aperture=centre, bridge=bridge_stats if len(centers)==2 else None,
                background_wavelength_test=background, background_suspect=abs(background["sn"] or 0)>=5,
                annular_psf_excess=bridge_stats if len(centers)==1 else None,
                seeing_arcsec=seeing, pixel_scale_arcsec=scale, aperture_radius_arcsec=1.5,
                ra=float(ra), dec=float(dec), two_psf_beta=beta,
                bridge_spaxels=int(bridge.sum()) if len(centers)==2 else 0,
                line_window_air_A=[float(wave[pi-half]),float(wave[pi+half])],
                mask="all nonzero flags excluded", aperture_unit="1e-17 erg/s/cm2",
                central_above_5sigma=(centre["sn"] or 0)>=5,
                annular_above_5sigma=len(centers)==1 and (bridge_stats["sn"] or 0)>=5,
                bridge_above_5sigma=len(centers)==2 and (bridge_stats["sn"] or 0)>=5)

    if halo_trials:
        result["halo_trials"] = []
        if centre["scatter"] and bridge_stats["scatter"] and len(centers)==1:
            from scipy.ndimage import gaussian_filter
            # Spectral injection is a sigma-2.4 A Gaussian truncated beyond
            # 4 sigma. None enters the sidebands or offset control windows.
            bins = np.arange(-4,5)*2.
            line = np.exp(-0.5*(bins/2.4)**2)
            fraction = np.sum(line[np.abs(bins)<=6])/np.sum(line)
            psf = psfs[0]/psfs[0].sum()
            for radius in (0.,1.5,3.,5.):
                halo = gaussian_filter(psf,radius/scale,mode="constant",cval=0.) if radius else psf
                # Total flux reference is the infinite-sky point/extended input.
                # Leakage of the convolved halo off the finite cube is retained.
                profile = 0.35*psf + 0.65*halo if radius else psf
                for nominal in (0.,5.,10.,20.,40.):
                    if not nominal and radius:
                        continue
                    total = nominal * centre["scatter"] / (fraction*np.sum(psf[ap]))
                    added = fraction*total*profile
                    cval, aval, _ = measure(m+added,good)
                    cstat, astat = empirical(cval,offset[:,0]),empirical(aval,offset[:,1])
                    result["halo_trials"].append(dict(
                        intrinsic_gaussian_sigma_arcsec=radius, halo_fraction=0.65 if radius else 0.0,
                        nominal_point_source_aperture_snr=nominal, total_flux=total*1e-17,
                        annular_surface_brightness=float(np.mean(added[region])*1e-17/scale**2),
                        central_snr=cstat["sn"], annular_snr=astat["sn"],
                        annular_passed=(astat["sn"] or 0)>=5,
                        central_and_annular_passed=(cstat["sn"] or 0)>=5 and (astat["sn"] or 0)>=5))
    return result


def blank_cube(file, positions, *, halo_meta=None, trial_output=None):
    """1.5-arcsec-radius apertures ≥5 arcsec from catalogued detections.

    Independent spatial placement, same-wavelength empirical error inflation,
    and conservative removal of every quality bit limit false precision from
    correlations introduced by spatial and wavelength interpolation.
    """
    with fits.open(file, memmap=False) as h:
        hd = h["DATA"]
        w = WCS(hd.header).celestial
        scale = float(np.mean(proj_plane_pixel_scales(w))*3600)
        xbad,ybad = w.world_to_pixel_values(positions[:,0],positions[:,1])
        yy,xx = np.mgrid[:hd.shape[1],:hd.shape[2]]
        flux,err = hd.data,h["ERROR"].data
        mask = h["MASK"].data if "MASK" in h else h["BITMASK"].data
        centers, fs, es = [],[],[]
        spacing = int(np.ceil(5.0/scale))
        for cy in range(spacing,hd.shape[1]-spacing,spacing):
            for cx in range(spacing,hd.shape[2]-spacing,spacing):
                if np.any(np.hypot(xbad-cx,ybad-cy)*scale<5.0):
                    continue
                ap = np.hypot(xx-cx,yy-cy)*scale<1.5
                a,e = flux[:,ap],err[:,ap]
                ok = np.isfinite(a) & (e>0) & np.isfinite(e) & (mask[:,ap]==0)
                good_wave = np.all(ok,axis=1)
                if np.mean(good_wave)<0.65:
                    continue
                # Cube unit is per 2-A bin. Transfer to per-A line-finder input.
                fs.append(np.where(good_wave, a.sum(axis=1)/2, np.nan))
                es.append(np.where(good_wave, np.sqrt((e**2).sum(axis=1))/2,np.nan))
                centers.append([float(v) for v in w.pixel_to_world_values(cx,cy)])
        if halo_meta is not None and trial_output is not None:
            trial_output["halo_injections"] = []
            midpoint = np.array(w.pixel_to_world_values((hd.shape[2]-1)/2,(hd.shape[1]-1)/2))
            # Select the blank nearest the IFU centre without inspecting
            # line flux or searching for a favorable wavelength or mask.
            if centers:
                xy = np.array(centers)
                dx = (xy[:,0]-midpoint[0]+180)%360-180
                i = int(np.argmin((dx*np.cos(np.deg2rad(midpoint[1])))**2+(xy[:,1]-midpoint[1])**2))
                for lam in (3950.,4850.):
                    rec = dict(ra=centers[i][0],dec=centers[i][1],wave=lam,
                               fwhm=halo_meta.get("fwhm"))
                    test = _cube_measure(h, rec, Path(file).name, halo_trials=True,resident=True)
                    trial_output["halo_injections"].append(dict(test,lambda_air=lam))
        if len(fs)<12:
            return [],np.empty((0,1036)),np.empty((0,1036))
        fs,es = np.array(fs),np.array(es)
        standard = fs/es
        med = np.nanmedian(standard,axis=0)
        inflation = 1.4826*np.nanmedian(np.abs(standard-med),axis=0)
        # Smooth the estimated spatial covariance in 22-A neighborhoods.
        from scipy.ndimage import median_filter
        inflation = median_filter(np.where(np.isfinite(inflation),inflation,1.),size=11)
        return centers, fs, es*np.maximum(inflation,1.)


def halo_summary(out):
    """Aggregate per-radius/flux-bin rates; trials sharing a blank are dependent."""
    import gzip
    import json
    from collections import defaultdict
    buckets = defaultdict(list)
    locs = 0
    visits = 0
    failed = []
    with gzip.open(Path(out)/"offsource.jsonl.gz","rt") as source:
        for entry in source:
            meta = json.loads(entry)["input"]
            visits += 1
            for loc in meta.get("halo_injections",[]):
                trials = loc.get("halo_trials",[])
                if not trials:
                    failed.append({"cube":loc.get("cube"),"wave":loc["lambda_air"],"status":loc["status"]})
                    continue
                locs += 1
                for trial in trials:
                    label = (trial["intrinsic_gaussian_sigma_arcsec"],trial["nominal_point_source_aperture_snr"])
                    buckets[label].append(trial)
    bins = []
    for (radius, nominal),trials in sorted(buckets.items()):
        bins.append({"intrinsic_sigma_arcsec":radius,"nominal_point_aperture_snr":nominal,
                     "n":len(trials),"annular_pass":sum(t["annular_passed"] for t in trials),
                     "central_and_annular_pass":sum(t["central_and_annular_passed"] for t in trials),
                     "minmax_total_flux":[float(fn([t["total_flux"] for t in trials])) for fn in (min,max)],
                     "minmax_annular_surface_brightness":[float(fn([t["annular_surface_brightness"] for t in trials]))
                                                          for fn in (min,max)]})
    (Path(out)/"halo-injections-summary.json").write_text(json.dumps(
        {"cube_attempts":visits,"usable_position_wavelength_pairs":locs,"unmeasured_pairs":failed,
         "bins":bins},indent=2)+"\n")
