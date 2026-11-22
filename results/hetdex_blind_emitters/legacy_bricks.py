"""Three-arcsecond positional checks against native LS DR9/DR10 Tractor bricks.

Layouts, flux and mask definitions: https://www.legacysurvey.org/dr10/files/
and https://www.legacysurvey.org/dr9/files/ . A zero-length catalog cone does
not constrain continuum brightness. The native r NEXP image is checked at the
actual sky coordinate, independently of catalog counterpart positions.
"""
import hashlib
import json
import math
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

ROOTS = ("https://portal.nersc.gov/cfs/cosmo/data/legacysurvey",
         "https://portal.nersc.gov/project/cosmo/data/legacysurvey")
CAPS = (("dr10", "south"), ("dr9", "south"), ("dr9", "north"))
RADIUS = 3.0


def _fetch(root, relpath, ceiling):
    """Preserve bytes and each unsuccessful attempt; never treat timeout as 404."""
    dest = root / relpath
    meta = dest.with_name(dest.name + ".receipt.json")
    if dest.exists() and meta.exists():
        receipt = json.loads(meta.read_text())
        if dest.stat().st_size == receipt["bytes"]:
            return dest, receipt
    dest.parent.mkdir(exist_ok=True, parents=True)
    attempts = []
    for host in ROOTS:
        url = host + "/" + relpath
        attempt = {"url": url}
        attempts.append(attempt)
        tmp = None
        try:
            with requests.get(url, stream=True, timeout=(5, 35),
                              headers={"Accept-Encoding": "identity"}) as response:
                attempt["http"] = response.status_code
                response.raise_for_status()
                if response.status_code != 200:
                    raise ValueError("incomplete HTTP response")
                length = response.headers.get("Content-Length")
                if length and int(length) > ceiling:
                    raise ValueError("object exceeds bounded download size")
                checksum = hashlib.sha256()
                count = 0
                with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as out:
                    tmp = Path(out.name)
                    for chunk in response.iter_content(256 * 1024):
                        count += len(chunk)
                        if count > ceiling:
                            raise ValueError("object exceeds bounded download size")
                        checksum.update(chunk)
                        out.write(chunk)
                if length and not response.headers.get("Content-Encoding") and count != int(length):
                    raise ValueError("truncated HTTP response")
            with tmp.open("rb") as stream:
                magic = stream.read(8)
            if magic != b"SIMPLE  " and not (relpath.endswith(".gz") and magic.startswith(b"\x1f\x8b")):
                raise ValueError("input lacks FITS or gzip signature")
            rec = dict(status="retrieved", path=relpath, url=url, bytes=count,
                       sha256=checksum.hexdigest(), attempts=attempts)
            os.replace(tmp, dest)
            meta.write_text(json.dumps(rec, allow_nan=False))
            return dest, rec
        except (ValueError, OSError, requests.RequestException) as exc:
            attempt["error"] = str(exc)
        finally:
            if tmp:
                tmp.unlink(missing_ok=True)
    missing = all(a.get("http") == 404 for a in attempts)
    return None, dict(path=relpath, status="not_published" if missing else "unavailable",
                      attempts=attempts)


def _geometry(file, points):
    rad = RADIUS / 3600.
    out = []
    with fits.open(file, memmap=False) as hdul:
        t = hdul[1].data
        if len(t) < 600000:
            raise ValueError("global brick grid is incomplete")
        name, lo, hi, bottom, top = (t[k] for k in ("brickname", "ra1", "ra2", "dec1", "dec2"))
        width = np.where(hi >= lo, hi - lo, 360 + hi - lo)
        for ra, dec in points:
            dx = np.mod(ra - lo, 360)
            pad = min(180., rad / max(1e-12, math.cos(math.radians(min(90., abs(dec) + rad)))))
            inside = (dx <= width) & (bottom <= dec) & (top >= dec)
            touch = ((dx <= width + pad) | (dx >= 360 - pad)) & (bottom <= dec + rad) & (top >= dec - rad)
            if not np.any(inside):
                raise ValueError("no sky brick contains this position")
            out.append([(_value(name[i]), bool(inside[i])) for i in np.flatnonzero(touch)])
    return out


def _value(value):
    value = value.item() if isinstance(value, np.generic) else value
    if isinstance(value, bytes):
        return value.decode("ascii").strip()
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _tractor(file, ra, dec):
    with fits.open(file, memmap=True) as hdul:
        t = hdul[1].data
        cols = {k.lower(): k for k in t.names}
        required = {"ra", "dec", "brick_primary", "brickid", "objid", "release",
                    "type", "maskbits", "flux_r", "flux_ivar_r", "nobs_r", "fracmasked_r", "fracflux_r"}
        if not required.issubset(cols):
            raise ValueError("missing native Tractor columns: " + str(sorted(required - cols.keys())))
        rar, decr = np.deg2rad(t[cols["ra"]]), np.deg2rad(t[cols["dec"]])
        hav = np.sin((decr - math.radians(dec))/2)**2 + np.cos(decr)*math.cos(math.radians(dec))*np.sin((rar-math.radians(ra))/2)**2
        sep = np.rad2deg(2 * np.arcsin(np.sqrt(np.clip(hav, 0., 1.)))) * 3600
        fields = set(required) | {"brickname", "fitbits", "ls_id", "wise_coadd_id"}
        for band in ("g", "r", "i", "z", "w1", "w2"):
            for prefix in ("flux", "flux_ivar", "fracflux", "fracmasked", "mw_transmission",
                           "nobs", "anymask", "allmask", "psfdepth"):
                fields.add(prefix + "_" + band)
        result = []
        # Each brick overlaps its neighbours; BRICK_PRIMARY assigns the unique owner.
        for i in np.flatnonzero((sep <= RADIUS) & np.asarray(t[cols["brick_primary"]], bool)):
            row = {k: _value(t[cols[k]][i]) for k in sorted(fields & cols.keys())}
            row.update(sep_arcsec=float(sep[i]),
                       tractor_id=f"{row['release']}-{row['brickid']}-{row['objid']}")
            for band in ("g", "r", "i", "z", "w1", "w2"):
                iv = row.get("flux_ivar_" + band)
                row["flux_err_" + band] = 1 / math.sqrt(iv) if iv and iv > 0 else None
            trans = row.get("mw_transmission_r")
            tests = {"unique_model": row["type"] != "DUP",
                     "maskbits_zero": row["maskbits"] == 0, "r_exposures": (row["nobs_r"] or 0) > 0,
                     "r_flux_and_ivar": row["flux_r"] is not None and row["flux_err_r"] is not None,
                     "low_masked_fraction": row["fracmasked_r"] is not None and 0 <= row["fracmasked_r"] < .15,
                     "low_contamination": row["fracflux_r"] is not None and 0 <= row["fracflux_r"] < .5,
                     "valid_extinction": trans is not None and 0 < trans <= 1.001,
                     "r_allmask_clear": row.get("allmask_r") == 0}
            row.update(quality_r=all(tests.values()), failed_r_quality=[k for k in tests if not tests[k]])
            result.append(row)
    return result


def _nexp(file, ra, dec):
    with fits.open(file, memmap=False) as hdul:
        for h in hdul:
            if h.header.get("NAXIS") != 2:
                continue
            w = WCS(h.header).celestial
            if w.pixel_n_dim != 2:
                continue
            xy = w.all_world2pix([[ra, dec]], 0)[0]
            if not np.isfinite(xy).all():
                continue
            x, y = np.floor(xy + .5).astype(int)
            if 0 <= x < h.header["NAXIS1"] and 0 <= y < h.header["NAXIS2"]:
                value = float(h.section[y:y+1, x:x+1][0, 0])
                if not math.isfinite(value) or value < 0 or int(value) != value:
                    raise ValueError("invalid number of exposures")
                return dict(status="measured", nexp_r=int(value), x0=int(x), y0=int(y), ra=ra, dec=dec)
    raise ValueError("position outside WCS grid")


def _paths(release, cap, brick):
    base = f"{release}/{cap}"
    coadd = f"{base}/coadd/{brick[:3]}/{brick}/"
    return (f"{base}/tractor/{brick[:3]}/tractor-{brick}.fits",
            coadd, coadd + f"legacysurvey-{brick}-nexp-r.fits.fz")


def native_cache_directory(cache):
    directory = cache if isinstance(cache, (str, os.PathLike)) else cache.root
    return Path(directory) / "legacy-tractor"


def legacy_brick_cones(points, cache):
    """Return one full evidence object per input coordinate, preserving order."""
    coords = [(float(p["ra"]), float(p["dec"])) if isinstance(p, dict)
              else (float(p[0]), float(p[1])) for p in points]
    if not coords:
        return []
    root = native_cache_directory(cache)
    geometry_receipts = []
    touched = None
    for release in ("dr10", "dr9"):
        p, receipt = _fetch(root, f"{release}/survey-bricks.fits.gz", 64 << 20)
        geometry_receipts.append(receipt)
        if p:
            try:
                touched = _geometry(p, coords)
                break
            except Exception as exc:
                receipt["geometry_error"] = str(exc)
    if touched is None:
        return [dict(status="unavailable", outcome="unavailable", ra=ra, dec=dec, matches=[],
                     radius_arcsec=RADIUS, selected_release=None, coverage=[], receipts=geometry_receipts,
                     r_exposure=dict(status="unavailable"), error="unverified sky brick grid") for ra, dec in coords]

    keys = {b for group in touched for b, _ in group}
    paths = [_paths(rel, cap, b)[0] for rel, cap in CAPS for b in sorted(keys)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        files = dict(zip(paths, pool.map(lambda s: _fetch(root, s, 192 << 20), paths)))
        exposures = {_paths(rel, cap, b)[2] for rel, cap in CAPS for group in touched for b, center in group
                     if center and files[_paths(rel, cap, b)[0]][0]}
        expfiles = dict(zip(sorted(exposures), pool.map(lambda s: _fetch(root, s, 48 << 20), sorted(exposures))))
    output = []
    for (ra, dec), group in zip(coords, touched, strict=True):
        receipts, coverage, candidates = {}, [], []
        latitude = SkyCoord(ra=ra*u.deg, dec=dec*u.deg).galactic.b.deg
        preferred = (("dr9","north"), *CAPS[:2]) if dec >= 32.375 and latitude > 0 else CAPS
        for release, cap in preferred:
            entries, matches, exp = [], [], []
            for b, center in group:
                path, coadd, ipath = _paths(release, cap, b)
                local, receipt = files[path]
                receipts[path] = receipt
                entry = dict(release=release, cap=cap, brickname=b, contains_center=center,
                             tractor_url=receipt.get("url", ROOTS[0]+"/"+path), coadd_root=ROOTS[0]+"/"+coadd,
                             status=receipt["status"])
                if local:
                    try:
                        rows = _tractor(local, ra, dec)
                        matches.extend(rows)
                        entry.update(status="queried", primary_matches=len(rows))
                    except Exception as exc:
                        entry.update(status="unavailable", error=str(exc))
                    if center:
                        fl, ir = expfiles[ipath]
                        receipts[ipath] = ir
                        image = dict(status="unavailable", input=ipath, url=ir.get("url", ROOTS[0]+"/"+ipath),
                                     sha256=ir.get("sha256"), brickname=b, release=release, cap=cap)
                        if fl:
                            try:
                                image.update(_nexp(fl, ra, dec))
                            except Exception as exc:
                                image["error"] = str(exc)
                        entry["r_exposure"] = image
                        exp.append(image)
                entries.append(entry)
            coverage.extend(entries)
            if (any(e["contains_center"] and e["status"] == "queried" for e in entries)
                    and all(e["status"] in ("queried","not_published") for e in entries)):
                # Prefer an exposure at the source itself, then the standard cap, independent of flux.
                image = max(exp, key=lambda d: d.get("nexp_r",-1), default={"status":"unavailable"})
                priority = 0 if image.get("nexp_r",0) > 0 else (2 if image.get("nexp_r") == 0 else 1)
                candidates.append((priority, release, cap, matches, image))
        all_missing = all(e["status"] == "not_published" for e in coverage)
        selected = min(candidates, key=lambda c:c[0]) if candidates else None
        if selected:
            _, rel, cap, matches, image = selected
            matches.sort(key=lambda m:(m["sep_arcsec"],m["tractor_id"]))
            status, outcome = "queried", "matched" if matches else "empty_catalog_cone"
        else:
            rel, cap, matches, image = None, None, [], {"status":"unavailable"}
            status = "outside_footprint" if all_missing else "unavailable"
            outcome = "no_published_brick" if all_missing else "unavailable"
        output.append(dict(status=status, outcome=outcome, ra=ra, dec=dec, radius_arcsec=RADIUS,
                           selected_release=rel, selected_cap=cap, table="native Legacy Surveys Tractor",
                           matches=matches, flux_units="observed nanomaggies (including Galactic attenuation)",
                           coverage=coverage, receipts=geometry_receipts+list(receipts.values()), r_exposure=image))
    return output
