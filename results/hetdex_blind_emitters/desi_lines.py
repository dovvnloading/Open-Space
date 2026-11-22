"""Independent optical spectral tests in the DESI follow-up VAC and NOIRLab SPARCL.

Flux fits below are integrated physical fluxes through the external instrument's fibre.
Each template is a Gaussian density with unit integral. Nearby components are solved
jointly with a linear local continuum. The secondary grid brackets modest centroid
offsets and linewidths, including the Ly-alpha to UV systemic velocity difference.
"""

import csv
from hashlib import sha256
import io
import json
from pathlib import Path

import numpy as np
import requests


C = 299792.458
TAP = "https://datalab.noirlab.edu/tap/sync"
SPARCL = "https://astrosparcl.datalab.noirlab.edu/api/spectras/"
VACUUM = {
    "lya": 1215.6701, "heii": 1640.42, "civ1548": 1548.204, "civ1551": 1550.781,
    "oiii1661": 1660.809, "oiii1666": 1666.150, "ciii1907": 1906.683,
    "ciii1909": 1908.734, "alt_oii3727": 3727.092, "alt_oii3730": 3729.875,
    "alt_neiii": 3869.860, "alt_hbeta": 4862.683, "alt_oiii4959": 4960.295,
    "alt_oiii5007": 5008.240, "alt_halpha": 6564.614,
}


def air_to_vac(wave):
    s2 = (1e4 / wave) ** 2
    return wave * (1.0 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


def fits_for_group(wave, flux, err, centers, velocities, widths):
    if min(centers) - 32 < wave[0] or max(centers) + 32 > wave[-1]:
        return "edge", []
    idx = (wave >= min(centers) - 45) & (wave <= max(centers) + 45)
    x, y, e = wave[idx], flux[idx], err[idx]
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0) & (y != -999)
    if valid.sum() < 20:
        return "masked", []
    # Masked cores, survey gaps and sky-line holes cannot produce upper limits.
    x, y, e = x[valid], y[valid], e[valid]
    possible = np.unique(wave[idx])
    dx = np.median(np.diff(possible))
    fits = []
    for dv in velocities:
        trial = np.asarray(centers) * (1 + dv / C)
        for width in widths:  # observed Gaussian sigma including instrumental broadening
            sigma = trial * width / C
            if any(np.sum(abs(x - c) < 1.7*s) < max(2, 0.60*3.4*s/dx)
                   for c, s in zip(trial, sigma)):
                continue
            basis = [np.ones_like(x), (x - np.mean(centers))/40]
            basis.extend(np.exp(-0.5*((x-c)/s)**2) / (np.sqrt(2*np.pi)*s)
                         for c, s in zip(trial, sigma))
            weighted = np.array(basis).T / e[:, None]
            normal = weighted.T @ weighted
            if np.linalg.cond(normal) > 1e11:
                continue
            cov = np.linalg.inv(normal)
            coef = cov @ (weighted.T @ (y/e))
            resid = y/e - weighted @ coef
            # Use the empirical local dispersion only to inflate nominal uncertainty.
            inflate = max(1., 1.4826*np.median(abs(resid - np.median(resid))))
            fits.append((dv, width, coef[2:] * 1e-17, cov[2:, 2:] * (inflate*1e-17)**2,
                         float(resid @ resid)))
    return ("ok" if fits else "masked"), fits


def measure_arrays(wave, flux, err, wave_air, channel):
    wave, flux, err = (np.asarray(a, dtype=float).ravel() for a in (wave, flux, err))
    if not (len(wave) == len(flux) == len(err) and np.all(np.diff(wave) > 0)):
        raise ValueError("External spectrum lengths or wavelength order mismatch")
    observed = float(air_to_vac(float(wave_air)))
    channel = str(channel).lower()
    primary = "heii" if "heii" in channel else "lya" if channel in ("ew", "blob", "lya", "lae") else None
    z = observed/VACUUM[primary] - 1 if primary else None
    z_oii = observed/3728.4835 - 1
    lines = {}
    def group(keys, centers, total=None, allow_offset=True):
        offsets = (-300, -150, 0, 150, 300) if allow_offset else (-100, 0, 100)
        status, fits = fits_for_group(wave, flux, err, centers, offsets, (65, 100, 175))
        nominal = next((fit for fit in fits if fit[:2] == (0, 100)), None)
        group_keys = list(keys) + ([total] if total else [])
        for j, label in enumerate(group_keys):
            result = {"wave_vac": float(centers[j] if j < len(centers) else np.mean(centers)),
                      "coverage": status if nominal is not None else ("edge" if status=="edge" else "masked")}
            if nominal:
                a = np.eye(len(keys))[j] if j < len(keys) else np.ones(len(keys))
                value = float(a @ nominal[2])
                sigma = float(np.sqrt(a @ nominal[3] @ a))
                limits = [max(0., float(a @ f[2])) + 3*np.sqrt(a @ f[3] @ a) for f in fits]
                result.update(flux=value, err=sigma, snr=value/sigma, upper3sigma=max(0., value)+3*sigma,
                              grid_trials=len(fits), grid_expected=len(offsets)*3,
                              # Missing trials invalidate a full-window nondetection claim.
                              max_upper3sigma=float(max(limits)) if len(fits)==len(offsets)*3 else None,
                              peak_snr=float(max(a @ f[2] / np.sqrt(a @ f[3] @ a) for f in fits)))
            lines[label] = result

    group(["primary"], [observed], allow_offset=False)
    if primary:
        for keys, total in ((["lya"], None), (["heii"], None),
                            (["civ1548","civ1551"], "civ"), (["oiii1661","oiii1666"], "oiii"),
                            (["ciii1907","ciii1909"], "ciii")):
            group(keys, [VACUUM[k]*(1+z) for k in keys], total)
    if z_oii >= 0:
        keys = ["alt_oii3727", "alt_oii3730"]
        group(keys, [VACUUM[k]*(1+z_oii) for k in keys], "alt_oii", allow_offset=False)
        for key in ("alt_neiii", "alt_hbeta", "alt_oiii4959", "alt_oiii5007", "alt_halpha"):
            group([key], [VACUUM[key]*(1+z_oii)], allow_offset=False)
    return {"assumed_primary":primary, "z_topocentric":z, "alternative_oii_z":z_oii,
            "flux_unit":"erg s-1 cm-2", "lines":lines,
            "method":"100 km/s observed sigma, fixed centers, local linear continuum; joint doublets, "
                     "summed covariance for blends. Grid covers 65,100,175 km/s with centroid shifts "
                     "up to 300 km/s for UV lines, 100 km/s for primary and low-z tests. Peak S/N "
                     "is descriptive across trials. Limits apply only to narrow emission in the "
                     "external fibre and require complete clean coverage of the template."}


def desi_lines(hdul, row, wave_air, channel):
    """Read DESI_HETDEX_SPEC_v1.6.fits while accepting either on-disk axis order."""
    nrow = len(hdul["INFO"].data)
    def vector(name):
        a = hdul[name].data
        if a.ndim == 1 or a.shape.count(1):
            return a.ravel()
        if a.shape[0] == nrow:
            return a[row]
        if a.shape[1] == nrow:
            return a[:, row]
        raise ValueError(f"Unknown {name} shape: {a.shape}")
    return measure_arrays(vector("DESI_WAVE_VACUUM"), vector("DESI_SPEC"),
                          vector("DESI_SPEC_ERR"), wave_air, channel)


def external_lines(ra, dec, wave_air, channel, cache):
    """A 1.25 arcsec spectroscopic crossmatch. Keep nonmatches distinct from outages."""
    ra, dec = float(ra), float(dec)
    radius = 1.25
    query = (f"SELECT sparcl_id,ra,dec,data_release,specid,spectype,redshift,redshift_warning "
             f"FROM sparcl.main WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
             f"CIRCLE('ICRS',{ra:.9f},{dec:.9f},{radius/3600:.10f}))")
    result = {"radius_arcsec":radius, "tap":TAP, "query":query, "matches":[]}
    cache = Path(cache)/"external-spectra"
    cache.mkdir(parents=True, exist_ok=True)
    try:
        cat = requests.get(TAP, params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":query},
                           timeout=(10,40))
        cat.raise_for_status()
        result["catalog_sha256"] = sha256(cat.content).hexdigest()
        table = csv.DictReader(io.StringIO(cat.text))
        if not table.fieldnames or "sparcl_id" not in [x.strip().lower() for x in table.fieldnames]:
            raise ValueError("Unexpected TAP response: " + cat.text[:240])
        for raw in table:
            obj = {k.lower():v for k,v in raw.items()}
            sep = float(2*np.arcsin(np.sqrt(np.sin(np.radians(float(obj["dec"])-dec)/2)**2 +
                        np.cos(np.radians(dec))*np.cos(np.radians(float(obj["dec"])))*
                        np.sin(np.radians(float(obj["ra"])-ra)/2)**2))*180/np.pi*3600)
            if sep > radius + 0.002:
                continue
            obj["sep_arcsec"] = sep
            uid = obj["sparcl_id"]
            if not (obj["data_release"].startswith(("DESI-", "SDSS-", "BOSS-"))):
                obj["status"] = "outside_datasets"
                result["matches"].append(obj)
                continue
            try:
                # API v14 supports a plain JSON response with metadata followed by records.
                params = {"format":"json","units":"False", "include":"wavelength,flux,ivar,mask",
                          "dataset_list":obj["data_release"]}
                local = cache / (sha256(uid.encode()).hexdigest() + ".json")
                if local.exists():
                    data = local.read_bytes()
                else:
                    response = requests.post(SPARCL, json=[uid], params=params, timeout=(10,45))
                    response.raise_for_status()
                    data = response.content
                    parsed = json.loads(data)
                    if not isinstance(parsed,list) or len(parsed) < 2:
                        raise ValueError("SPARCL lookup returned no spectra")
                    local.write_bytes(data)
                record = json.loads(data)[1]
                w, f, iv = (np.asarray(record[k],float) for k in ("wavelength","flux","ivar"))
                with np.errstate(invalid="ignore", divide="ignore"):
                    e = np.where(iv > 0, 1/np.sqrt(iv), np.nan)
                if "mask" in record:
                    e[np.asarray(record["mask"])!=0] = np.nan
                obj.update(status="measured", spectra_url=SPARCL, spectral_sha256=sha256(data).hexdigest(),
                           diagnostics=measure_arrays(w, f, e, wave_air, channel))
            except (OSError, ValueError, KeyError, IndexError, TypeError) as err:
                obj["status"] = "spectrum_unavailable"
                obj["reason"] = str(err)[:320]
            result["matches"].append(obj)
        result["status"] = "matched" if result["matches"] else "no_counterpart"
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as err:
        result.update(status="unavailable", reason=str(err)[:320])
    return result
