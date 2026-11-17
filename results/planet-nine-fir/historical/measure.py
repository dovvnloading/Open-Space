"""Reproduce historical catalogue checks and coadd photometry. Python, numpy, astropy."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
IRAS = SkyCoord(35.74075, -48.5125, unit="deg")
AKARI = SkyCoord(35.18379, -49.2135, unit="deg")
CONTROL = SkyCoord(34.05346, -49.45566, unit="deg")  # AKARI BSC J0216128-492720
GALAXY = SkyCoord(35.7544007, -48.5112676, unit="deg")
TARGETS = {"IRAS endpoint": IRAS, "AKARI endpoint": AKARI,
           "stationary control": CONTROL, "nearby W4 galaxy": GALAXY}


def load_inputs(cache):
    inputs = json.loads((ROOT / "inputs.json").read_text())
    cache.mkdir(parents=True, exist_ok=True)
    for item in inputs:
        path = cache / item["name"]
        if not path.exists():
            with urlopen(item["url"], timeout=120) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise ValueError(f"Archive checksum changed for {path.name}")
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"Checksum mismatch for {path}")
    return inputs


def catalogue_checks(cache):
    with gzip.open(cache / "cfa_cds_fsc", "rt") as file:
        row = next(r for r in file if r.startswith("F02211-4844"))
    ra = 15 * (int(row[12:14]) + int(row[14:16])/60 + int(row[16:19])/36000)
    dec = (int(row[20:22]) + int(row[22:24])/60 + int(row[24:26])/3600) * (-1 if row[19] == "-" else 1)
    fsc_position = SkyCoord(ra, dec, unit="deg", frame="fk4", equinox="B1950", obstime="J1983.5").icrs
    # Restrict by declination before computing exact great-circle distances.
    with gzip.open(cache / "bsc", "rt") as file:
        nearby = [r for r in file if -53 < float(r[36:45]) < -45 and 28 < float(r[26:35]) < 43]
    coords = SkyCoord([float(r[26:35]) for r in nearby], [float(r[36:45]) for r in nearby], unit="deg")
    nearest = {}
    for name, coord in TARGETS.items():
        distances = coords.separation(coord).arcsec
        idx = int(np.argmin(distances))
        nearest[name] = {"name": nearby[idx][11:25], "separation_arcsec": float(distances[idx]),
                         "within_60arcsec": int(np.count_nonzero(distances <= 60))}
    return {"fsc_line": row.rstrip(), "fsc_icrs_deg": [fsc_position.ra.deg, fsc_position.dec.deg],
            "fsc_minus_paper_separation_arcsec": fsc_position.separation(IRAS).arcsec,
            "two_paper_positions_separation_arcmin": IRAS.separation(AKARI).arcmin,
            "FSC_60_100_Jy": [float(row[65:74]), float(row[74:83])],
            "FSC_60_100_quality": [int(row[85]),int(row[86])],
            "FSC_60_100_error_pct": [int(row[93:96]), int(row[96:99])],
            "FSC_major_minor_arcsec_PA_deg": [int(row[26:29]),int(row[29:32]),int(row[32:35])],
            "BSCv1_nearest": nearest}


def pixel_position(wcs, sky, akari):
    spherical = sky.barycentricmeanecliptic.spherical if akari else sky.fk4.spherical
    return tuple(float(v) for v in wcs.all_world2pix(spherical.lon.deg, spherical.lat.deg, 0))


def forced(data, xy, pixel_arcsec, fwhm, aperture_correction=None):
    """Gaussian+plane for IRAS, or the Takita et al. aperture for AKARI."""
    x, y = xy
    radius = 300 if aperture_correction is not None else 2*fwhm
    n = int(np.ceil(radius/pixel_arcsec)) + 1
    ix, iy = round(x), round(y)
    if min(ix,iy) < n or ix+n >= data.shape[1] or iy+n >= data.shape[0]:
        return np.nan
    cut = data[iy-n:iy+n+1,ix-n:ix+n+1]
    dx = (np.arange(ix-n,ix+n+1)-x)*pixel_arcsec
    dy = (np.arange(iy-n,iy+n+1)-y)*pixel_arcsec
    a,b = np.meshgrid(dx,dy)
    r = np.hypot(a,b)
    good = np.isfinite(cut) & (cut > -1e4)
    if aperture_correction is not None:
        aperture = r <= 90
        annulus = (r >= 120) & (r <= 300)
        if not good[aperture].all() or good[annulus].mean() < .9:
            return np.nan
        return (cut[aperture]-np.median(cut[annulus & good])).sum() * (pixel_arcsec/206264.806)**2 * 1e6 / aperture_correction
    fit = r <= radius
    if good[fit].mean() < .9 or not good[r < fwhm/2].all():
        return np.nan
    fit &= good
    template = np.exp(-4*np.log(2)*(r[fit]/fwhm)**2)
    design = np.stack([template, np.ones_like(template),a[fit]/fwhm,b[fit]/fwhm],axis=1)
    amplitude = np.linalg.lstsq(design,cut[fit],rcond=None)[0][0]
    beam_sr = np.pi/(4*np.log(2))*(fwhm/206264.806)**2
    return amplitude * beam_sr * 1e6


def map_checks(label, data, header, band, coverage=None):
    data = np.squeeze(data)
    if data.ndim != 2:
        raise ValueError(f"Unexpected image shape {data.shape}")
    akari = label.startswith("AKARI")
    fwhm = {60:240,100:258,65:63.4,90:77.8}[band]
    corr = {65:.627,90:.696}.get(band)
    wcs = WCS(header, naxis=2)
    pixelscale = 3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
    out = []
    for name, coord in TARGETS.items():
        xy = pixel_position(wcs,coord,akari)
        flux = forced(data,xy,pixelscale,fwhm,corr)
        # Fixed grid of controls, approximately 10 arcmin spacing, 12--45 arcmin radii.
        # The same tangent-plane pixel scale is used for target and neighboring controls.
        controls = []
        for east in range(-40,41,10):
            for north in range(-40,41,10):
                if 12 <= np.hypot(east,north) <= 45:
                    pos = (xy[0]+east*60/pixelscale,xy[1]+north*60/pixelscale)
                    controls.append(forced(data,pos,pixelscale,fwhm,corr))
        valid = np.asarray(controls)
        valid = valid[np.isfinite(valid)]
        if len(valid) >= 5:
            med = np.median(valid)
            scatter = 1.4826*np.median(abs(valid-med))
        else:
            med,scatter = np.nan,np.nan
        ix,iy = map(round, xy)
        nscan = coverage[iy,ix] if coverage is not None and 0<=ix<data.shape[1] and 0<=iy<data.shape[0] else None
        row = {"product":label, "um":band, "position":name, "flux_Jy":flux, "blank_median_Jy":med,
               "blank_MAD_sigma_Jy":scatter, "contrast_MAD":(flux-med)/scatter if scatter>0 else np.nan,
               "n_blank":len(valid), "n_blank_at_or_above":int(np.count_nonzero(valid>=flux)) if np.isfinite(flux) else None,
               "center_nscan":nscan}
        out.append({k:(round(float(v),5) if np.isfinite(v) else None) if isinstance(v,(float,np.floating)) else v for k,v in row.items()})
    return out


def tar_fits(archive, member):
    blob = archive.extractfile(member).read()
    if member.endswith(".gz"):
        blob = gzip.decompress(blob)
    return fits.getdata(io.BytesIO(blob),header=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--output",type=Path,default=ROOT)
    args = ap.parse_args()
    load_inputs(args.cache)
    args.output.mkdir(parents=True,exist_ok=True)
    checks = catalogue_checks(args.cache)
    (args.output/"catalogues.json").write_text(json.dumps(checks, indent=2)+"\n")
    results = []
    for b,um in [(3,60),(4,100)]:
        for h in range(4):
            data,header = fits.getdata(args.cache/f"I045B{b}H{h}.fits",header=True)
            results.extend(map_checks(f"IRIS H{h}",data,header,um))
    for band,um in [("WideS",90),("N60",65)]:
        fn = f"l004.36_b-55.00_ecl_6deg_{band}"
        with tarfile.open(args.cache/(fn+".tar")) as archive:
            data,header = tar_fits(archive,fn+"_fixstripe.fits.gz")
            coverage,_ = tar_fits(archive,fn+"_nscan.fits.gz")
        results.extend(map_checks("AKARI AFASS",data,header,um,coverage))
    with tarfile.open(args.cache/"045.tar") as archive:
        for member in archive.getnames():
            name = Path(member).name
            if name.startswith("I045B") and name.endswith(".FIT"):
                band, hcon = int(name[5]), int(name[7])
                if band in (3,4):
                    data,header = tar_fits(archive,member)
                    results.extend(map_checks(f"ISSA H{hcon}",data,header,{3:60,4:100}[band]))
    with (args.output/"map_photometry.csv").open("w",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=results[0])
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
