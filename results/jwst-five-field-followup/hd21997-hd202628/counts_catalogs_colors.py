"""Catalog and ALMA associations, field counts and atmospheric grid at all blind triggers."""

import argparse
import csv
from pathlib import Path
import sys
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u

from auxiliary import FIELDS, GAIA, apparent


def rows(f):
    return list(csv.DictReader(f.open()))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("merged", type=Path)
    p.add_argument("cache", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    sys.path[:0] = [str(a.merged / "nircam"), str(a.merged / "miri")]
    from nircam_reduce import write_table
    from vet import aperture, models
    import vet

    bands = ("F200W", "F444W", "F1140C")
    vet.BANDS = bands
    pars, grid = models(a.cache / "bobcat.tar.gz")
    cats, counts, phot, col = [], [], [], []
    for field, group in FIELDS.items():
        sc = [s for s in rows(a.out / "auxiliary_peaks.csv") if s["field"] == field]
        for i, r in enumerate(rows(a.out / (field + "_peaks.csv")), 1):
            sc.append(
                dict(
                    id=i,
                    band="F444W",
                    file=field + "_F444W",
                    east=r["east_arcsec"],
                    north=r["north_arcsec"],
                    mjd=60463.09718005787,
                )
            )
        gcat = rows(here / "catalogs" / (field + "_GaiaDR3.csv"))
        host = next(r for r in gcat if r["source_id"] == GAIA[field])
        g0 = SkyCoord(float(host["ra"]), float(host["dec"]), unit="deg")
        catalogs = [
            ("GaiaDR3", gcat, "source_id", 2016.0),
            (
                "SIMBAD",
                rows(here / "catalogs" / (field + "_SIMBAD.csv")),
                "main_id",
                2000.0,
            ),
            (
                "2MASS",
                rows(here / "catalogs" / (field + "_2MASS.csv")),
                "designation",
                None,
            ),
        ]
        for label, table, key, epoch in catalogs:
            others = [
                r
                for r in table
                if not (key == "source_id" and r[key] == GAIA[field])
                and not ("HD" in r[key] and r[key].replace(" ", "") == field)
                and r[key] not in ("03315364-2536509", "21182725-4320048")
            ]
            for s in sc:
                m = float(s["mjd"])
                # Gaia catalog position is barycentric at J2016, so add the
                # apparent motion at each epoch including annual parallax.
                en = apparent(host, m) + np.asarray([s["east"], s["north"]], float)
                sky = g0.spherical_offsets_by(*(en * u.arcsec))
                dd = []
                for obj in others:
                    loc = SkyCoord(float(obj["ra"]), float(obj["dec"]), unit="deg")
                    if epoch and obj.get("pmra") and obj.get("pmdec"):
                        if label == "GaiaDR3" and obj.get("parallax"):
                            offset = apparent(obj, m)
                        else:
                            offset = (
                                np.asarray([obj["pmra"], obj["pmdec"]], float)
                                * (Time(m, format="mjd").jyear - epoch)
                                / 1000
                            )
                        loc = loc.spherical_offsets_by(*(offset * u.arcsec))
                    dd.append(sky.separation(loc).arcsec)
                i = int(np.argmin(dd)) if dd else None
                cats.append(
                    dict(
                        field=field,
                        file=s["file"],
                        id=s["id"],
                        band=s["band"],
                        ra=sky.ra.deg,
                        dec=sky.dec.deg,
                        catalog=label,
                        nearest_id=others[i][key] if dd else "",
                        separation_arcsec=dd[i] if dd else "",
                        inside_1arcsec=int(bool(dd) and dd[i] < 1),
                    )
                )
            # Cone catalogue counts explicitly exclude host. The control annulus
            # is fully within the 45" queried cone, including epoch shifts.
            loc = [
                SkyCoord(float(o["ra"]), float(o["dec"]), unit="deg") for o in others
            ]
            sep = np.array([g0.separation(z).arcsec for z in loc])
            nn = int(((sep >= 10) & (sep < 35)).sum())
            counts.append(
                dict(
                    field=field,
                    catalog=label,
                    control_inner_arcsec=10,
                    control_outer_arcsec=35,
                    control_count=nn,
                    scaled_to_2p5_5p5=nn * (5.5**2 - 2.5**2) / (35**2 - 10**2),
                )
            )
        files = [
            next(a.cache.glob("*" + group + "*" + b.lower() + "*i2d.fits"))
            for b in bands[:2]
        ]
        files += sorted(
            a.cache.glob(
                "*" + ("c1018" if field == "HD21997" else "c1017") + "_psfsub.fits"
            )
        )
        for f in files:
            with fits.open(f) as h:
                primary = h[0].header
                hdr = h["SCI"].header
                im = h["SCI"].data.astype(float) * hdr["PIXAR_SR"] * 1e12
                w = WCS(hdr, naxis=2)
                band = primary["FILTER"]
                mjd = primary["EXPSTART"]
                if im.ndim == 3:
                    from nircam_reduce import median

                    im[(h["DQ"].data & 3) != 0] = np.nan
                    im = median(im, axis=0)
                    cx, cy = hdr["CRPIX1"] - 1, hdr["CRPIX2"] - 1
                else:
                    im[h["WHT"].data <= 0] = np.nan
                    cx, cy = w.world_to_pixel_values(
                        primary["TARG_RA"], primary["TARG_DEC"]
                    )
                hostpix = w.pixel_to_world(cx, cy)
                rad = {"F200W": 2.5, "F444W": 2.5, "F1140C": 3.2}[band]
                cal = dict(radius=rad, skyin=rad + 2, skyout=rad + 5, apcorr=1 / 0.6)

                def at(xx, yy):
                    if (
                        min(xx, yy) < 12
                        or xx >= im.shape[1] - 12
                        or yy >= im.shape[0] - 12
                    ):
                        return np.nan
                    return aperture(im, xx, yy, cal)

                for s in sc:
                    key = dict(field=field, file=s["file"], id=s["id"])
                    en = np.array([s["east"], s["north"]], float)
                    difference = apparent(host, float(s["mjd"])) - apparent(host, mjd)
                    for track, offset in [
                        ("comoving", en),
                        ("stationary", en + difference),
                    ]:
                        tgt = hostpix.spherical_offsets_by(*(offset * u.arcsec))
                        xx, yy = w.world_to_pixel(tgt)
                        xx, yy = float(xx), float(yy)
                        phase = np.arange(2, 59) * 2 * np.pi / 60
                        controls = np.asarray(
                            [
                                at(
                                    cx + (xx - cx) * np.cos(t) - (yy - cy) * np.sin(t),
                                    cy + (xx - cx) * np.sin(t) + (yy - cy) * np.cos(t),
                                )
                                for t in phase
                            ]
                        )
                        sig = 1.4826 * np.nanmedian(
                            abs(controls - np.nanmedian(controls))
                        )
                        valid = (
                            abs(xx - cx) * np.sqrt(hdr["PIXAR_SR"]) * 206264.806 > 0.37
                            and abs(yy - cy) * np.sqrt(hdr["PIXAR_SR"]) * 206264.806
                            > 0.37
                            if band == "F1140C"
                            else True
                        )
                        flux = at(xx, yy)
                        phot.append(
                            dict(
                                **key,
                                sampled_file=f.name,
                                track=track,
                                x=xx,
                                y=yy,
                                band=band,
                                flux_proxy_uJy=flux,
                                control_MAD_uJy=sig,
                                usable=int(
                                    valid
                                    and np.isfinite(sig)
                                    and sig > 0
                                    and np.isfinite(flux)
                                ),
                            )
                        )
        for s in sc:
            same = [
                r
                for r in phot
                if r["field"] == field
                and r["file"] == s["file"]
                and r["id"] == s["id"]
                and r["track"] == "comoving"
                and r["usable"]
            ]
            flux, err, idx = [], [], []
            for j, b in enumerate(bands):
                rr = [r for r in same if r["band"] == b]
                if rr:
                    flux.append(np.mean([r["flux_proxy_uJy"] for r in rr]))
                    err.append(max(r["control_MAD_uJy"] for r in rr))
                    idx.append(j)
            result = dict(
                field=field,
                file=s["file"],
                id=s["id"],
                n_bands=len(idx),
                bands=",".join(bands[j] for j in idx),
                min_chisq="",
                teff_K="",
                logg="",
                scale="",
            )
            if len(idx) >= 2:
                v = np.array(flux)
                sigma = np.hypot(err, 0.2 * np.abs(v))
                g = grid[:, idx]
                fac = np.maximum(
                    0,
                    np.sum(g * v / sigma**2, axis=1) / np.sum((g / sigma) ** 2, axis=1),
                )
                chisq = np.sum(((fac[:, None] * g - v) / sigma) ** 2, axis=1)
                i = np.argmin(chisq)
                result.update(
                    min_chisq=chisq[i],
                    teff_K=pars[i, 0] if fac[i] > 0 else "",
                    logg=pars[i, 1] if fac[i] > 0 else "",
                    scale=fac[i],
                )
            col.append(result)
    for name, data in [
        ("catalog_matches", cats),
        ("catalog_background", counts),
        ("all_forced_colors", phot),
        ("all_model_comparisons", col),
    ]:
        write_table(a.out / (name + ".csv"), list(data[0]), data)


if __name__ == "__main__":
    main()
