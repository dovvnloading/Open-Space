"""Blind, per-visit MIRI/F200W checks, matched detections and forced-position measurements."""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
import astropy.units as u
from scipy import ndimage

FIELDS = {"HD21997": "c1007", "HD202628": "c1009"}
GAIA = {"HD21997": "5084963312478488960", "HD202628": "6577318075351824512"}
RA0 = {"HD21997": 52.97393, "HD202628": 319.61588}
DEC0 = {"HD21997": -25.61425, "HD202628": -43.33451}


def sky_offset(w, cx, cy, x, y):
    start = w.pixel_to_world(cx, cy)
    de, dn = start.spherical_offsets_to(w.pixel_to_world(x, y))
    return np.array([de.arcsec, dn.arcsec])


def apparent(gaia, mjd):
    ra, dec = np.radians([float(gaia["ra"]), float(gaia["dec"])])
    obs = get_body_barycentric("earth", Time(mjd, format="mjd")).xyz.to_value(u.au)
    proj = np.array(
        [
            [-np.sin(ra), np.cos(ra), 0],
            [-np.cos(ra) * np.sin(dec), -np.sin(ra) * np.sin(dec), np.cos(dec)],
        ]
    )
    motion = np.array([float(gaia["pmra"]), float(gaia["pmdec"])]) * (
        Time(mjd, format="mjd").jyear - 2016
    )
    return (motion - float(gaia["parallax"]) * (proj @ obs)) / 1000


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    sys.path[:0] = [str(a.merged / "miri"), str(a.merged / "nircam")]
    from search import smooth_response, peaks
    from nircam_reduce import median, write_table
    from vet import aperture, models
    import vet

    vet.BANDS = ("F200W", "F444W", "F1140C")
    model_par, model_flux = models(a.cache / "bobcat.tar.gz")
    det, vis, assoc, forced, colors = [], [], [], [], []
    for field, group in FIELDS.items():
        gaia = next(
            r
            for r in csv.DictReader(
                (here / "catalogs" / (field + "_GaiaDR3.csv")).open()
            )
            if r["source_id"] == GAIA[field]
        )
        base = json.loads((a.out / (field + "_measurement.json")).read_text())
        initial = list(csv.DictReader((a.out / (field + "_peaks.csv")).open()))
        if field == "HD21997":
            files = sorted(a.cache.glob("*c1018_psfsub.fits"))
        else:
            files = sorted(a.cache.glob("*c1017_psfsub.fits"))
        files += list(a.cache.glob("*" + group + "*f200w*i2d.fits"))
        by_file = {}
        for path in files:
            with fits.open(path) as hd:
                head = hd[0].header
                h = hd["SCI"].header
                w = WCS(h, naxis=2)
                im = hd["SCI"].data.astype(float)
                band = head["FILTER"]
                mjd = head["EXPSTART"]
                if im.ndim == 3:
                    im[(hd["DQ"].data & 3) != 0] = np.nan
                    im = median(im, axis=0)
                    cx, cy = h["CRPIX1"] - 1, h["CRPIX2"] - 1
                else:
                    im[hd["WHT"].data <= 0] = np.nan
                    cx, cy = w.world_to_pixel_values(head["TARG_RA"], head["TARG_DEC"])
                scale = np.sqrt(h["PIXAR_SR"]) * 206264.806
                y, x = np.indices(im.shape)
                r = np.hypot(x - cx, y - cy) * scale
                valid = np.isfinite(im)
                good = ndimage.minimum_filter(valid, size=5, mode="constant", cval=0)
                # Mask the two detector-fixed 4QPM boundaries by one FWHM on each side.
                if band == "F1140C":
                    good &= (abs(x - cx) * scale > 0.37) & (abs(y - cy) * scale > 0.37)
                search = good & (r >= 2.5) & (r < 5.5)
                filt = smooth_response(np.where(valid, im, np.nanmedian(im)))
                score = np.zeros_like(im)
                rings = np.arange(2.5, 5.51, 0.5)
                for lo, hi in zip(rings, rings[1:]):
                    ring = search & (r >= lo) & (r < hi)
                    if ring.sum() < 40:
                        continue
                    sample = filt[ring]
                    med = np.median(sample)
                    sig = 1.4826 * np.median(abs(sample - med))
                    score[ring] = (sample - med) / sig
                ps, sn = peaks(score, search)
                nneg = len(peaks(-score, search)[0])
                pos = [sky_offset(w, cx, cy, *pt) for pt in ps]
                this = []
                for k, (pt, en, s) in enumerate(zip(ps, pos, sn), 1):
                    row = dict(
                        field=field,
                        file=path.name,
                        id=k,
                        band=band,
                        mjd=mjd,
                        x=pt[0],
                        y=pt[1],
                        east=en[0],
                        north=en[1],
                        filter_snr=s,
                    )
                    det.append(row)
                    this.append(row)
                vis.append(
                    dict(
                        field=field,
                        file=path.name,
                        band=band,
                        mjd=mjd,
                        search_area_arcsec2=search.sum() * scale**2,
                        positive=len(pos),
                        negative=nneg,
                        masked_host_x=float(cx),
                        masked_host_y=float(cy),
                    )
                )
                by_file[path.name] = (mjd, this)
                for k, row in enumerate(initial, 1):
                    en = np.array(
                        [float(row["east_arcsec"]), float(row["north_arcsec"])]
                    )
                    for track in ("comoving", "stationary"):
                        en1 = (
                            en
                            if track == "comoving"
                            else en
                            + apparent(gaia, base["expstart_mjd"])
                            - apparent(gaia, mjd)
                        )
                        target = w.pixel_to_world(cx, cy).spherical_offsets_by(
                            *(en1 * u.arcsec)
                        )
                        xx, yy = w.world_to_pixel(target)
                        xx, yy = float(xx), float(yy)
                        ix, iy = round(xx), round(yy)
                        # Signed small-aperture flux proxy with local annular subtraction;
                        # model comparison uses 60% encircled energy estimate for each band.
                        rad = 2.5 if band == "F200W" else 3.2
                        cal = dict(
                            radius=rad, skyin=rad + 2, skyout=rad + 5, apcorr=1 / 0.6
                        )
                        im_ujy = im * h["PIXAR_SR"] * 1e12

                        def at(x, y):
                            if (
                                x < 12
                                or y < 12
                                or x > im.shape[1] - 12
                                or y > im.shape[0] - 12
                            ):
                                return np.nan
                            return aperture(im_ujy, x, y, cal)

                        phase = np.arange(60) * np.pi / 30
                        controls = np.array(
                            [
                                at(
                                    cx + (xx - cx) * np.cos(t) - (yy - cy) * np.sin(t),
                                    cy + (xx - cx) * np.sin(t) + (yy - cy) * np.cos(t),
                                )
                                for t in phase[2:-1]
                            ]
                        )
                        scatter = 1.4826 * np.nanmedian(
                            abs(controls - np.nanmedian(controls))
                        )
                        nearest = min(
                            (np.linalg.norm(np.asarray(p) - en1) for p in pos),
                            default=np.nan,
                        )
                        forced.append(
                            dict(
                                field=field,
                                peak=k,
                                track=track,
                                file=path.name,
                                mjd=mjd,
                                band=band,
                                x=xx,
                                y=yy,
                                east=en1[0],
                                north=en1[1],
                                valid=int(search[iy, ix]),
                                filter_snr=score[iy, ix],
                                nearest_6sigma_arcsec=nearest,
                                flux_proxy_uJy=at(xx, yy),
                                control_MAD_uJy=scatter,
                            )
                        )
        # Include every auxiliary detection, in each direction. A 0.35 arcsec
        # search radius admits registration/mask errors but may leave both
        # tracks unresolved. No orbit or significance inferred from these.
        by_file[field + "_F444W"] = (
            base["expstart_mjd"],
            [
                dict(id=i, east=float(r["east_arcsec"]), north=float(r["north_arcsec"]))
                for i, r in enumerate(initial, 1)
            ],
        )
        for f, (m, rows) in by_file.items():
            for g, (n, others) in by_file.items():
                if f == g:
                    continue
                d = apparent(gaia, m) - apparent(gaia, n)
                for s in rows:
                    en = np.array([s["east"], s["north"]])

                    def nearest(pred):
                        dist = [
                            np.linalg.norm(np.array([z["east"], z["north"]]) - pred)
                            for z in others
                        ]
                        return (
                            (others[np.argmin(dist)]["id"], min(dist))
                            if dist
                            else ("", np.nan)
                        )

                    i, dc = nearest(en)
                    j, ds = nearest(en + d)
                    assoc.append(
                        dict(
                            field=field,
                            from_file=f,
                            peak=s["id"],
                            to_file=g,
                            baseline_days=n - m,
                            hypothesis_separation_arcsec=np.linalg.norm(d),
                            comoving_nearest=i,
                            comoving_arcsec=dc,
                            stationary_nearest=j,
                            stationary_arcsec=ds,
                        )
                    )
        for k, row in enumerate(initial, 1):
            subset = [
                r
                for r in forced
                if r["peak"] == k
                and r["field"] == field
                and r["track"] == "comoving"
                and r["valid"]
            ]
            flux = [float(row["flux_uJy"])]
            err = [flux[0] / float(row["snr"])]
            bands = [1]
            for j, b in [(0, "F200W"), (2, "F1140C")]:
                r = [
                    x
                    for x in subset
                    if x["band"] == b
                    and np.isfinite(x["flux_proxy_uJy"])
                    and np.isfinite(x["control_MAD_uJy"])
                    and x["control_MAD_uJy"] > 0
                ]
                if r:
                    bands.append(j)
                    flux.append(np.mean([x["flux_proxy_uJy"] for x in r]))
                    err.append(max(x["control_MAD_uJy"] for x in r))
            flux, err = np.array(flux), np.array(err)
            err = np.hypot(err, 0.2 * np.abs(flux))
            grid = model_flux[:, bands]
            fac = np.maximum(
                0,
                np.sum(grid * flux / (err**2), axis=1)
                / np.sum((grid / err) ** 2, axis=1),
            )
            chi = np.sum(((grid * fac[:, None] - flux) / err) ** 2, axis=1)
            best = np.argmin(chi)
            colors.append(
                dict(
                    field=field,
                    peak=k,
                    bands=",".join(vet.BANDS[j] for j in bands),
                    flux_uJy=",".join(map(str, flux)),
                    effective_sigma_uJy=",".join(map(str, err)),
                    model_teff_K=model_par[best, 0],
                    logg=model_par[best, 1],
                    min_chisq=chi[best],
                    n_grid=len(grid),
                )
            )
    for name, rows in [
        ("auxiliary_peaks", det),
        ("auxiliary_visits", vis),
        ("motion_nearest", assoc),
        ("forced_tracks", forced),
        ("F444W_trigger_colors", colors),
    ]:
        write_table(a.out / (name + ".csv"), list(rows[0]), rows)


if __name__ == "__main__":
    main()
