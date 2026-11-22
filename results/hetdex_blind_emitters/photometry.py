"""Rest-frame Ly-alpha EW bounds from calibrated red-side imaging and HETDEX flux.

All quoted bounds condition on a power-law UV continuum with -3 <= beta <= 0.
Matched-aperture and model-total measurements remain separately identified.
"""

from functools import lru_cache
from hashlib import sha256
import math

import numpy as np
from scipy.integrate import trapezoid
from speclite.filters import load_filter


HSC_R = ("hsc2017-r", "hsc2017-r2")
LEGACY_R = ("BASS-r", "decam2014-r", "decamDR1-r")
NMAGGY_FNU = 3.631e-29  # 3631 Jy * 1e-9 * 1e-23 (erg/s/cm2/Hz)
C_A_PER_S = 2.99792458e18


def num(value, fallback=None):
    try:
        value = float(value)
        if np.isfinite(value) and value != -999:
            return value
    except (TypeError, ValueError):
        pass
    return fallback


@lru_cache(maxsize=12)
def response(name):
    f = load_filter(name)
    w, r = np.array(f.wavelength, dtype=float), np.array(f.response, dtype=float)
    digest = sha256(np.column_stack((w,r)).astype("<f8").tobytes()).hexdigest()
    return w, r, digest


def continuum_upper(lam, upper_nanomaggy, names):
    """Convert AB count rate to maximal f_lambda at lam, which is blue of r."""
    tests = []
    for name in names:
        w, r, digest = response(name)
        weight = r / w  # Photon-counting AB average weights f_nu by T(lambda)/lambda.
        full = trapezoid(weight, w)
        cut = lam + 25  # exclude Ly-alpha and any blueward forest absorption
        red = w > cut
        if red.sum() < 3:
            return None, tests
        red_fraction = float(trapezoid(weight[red], w[red])/full)
        tests.append({"name":name, "redward_weight":red_fraction, "curve_f64le_sha256":digest})
        if red_fraction < .99:
            return None, tests
        # Treat the entire excluded blue part as zero signal. For fixed f_lambda(lam)
        # and red wavelengths, the denominator is minimized at beta = -3.
        norm = trapezoid((w[red]/lam)**-3 * w[red] * r[red], w[red])
        tests[-1]["continuum_upper"] = float(upper_nanomaggy * NMAGGY_FNU * C_A_PER_S * full / norm)
    return max(row["continuum_upper"] for row in tests), tests


def bound(rec, origin, filt, photo_upper, line, line_err, aperture, eligible=True):
    z, lam = num(rec.get("z_lya")), num(rec.get("wave_vac"), num(rec.get("wave")))
    if z is None or lam is None or photo_upper <= 0 or line is None or not line_err or line_err <= 0:
        return None
    # Small absolute calibration margins are additional to quoted statistical errors.
    # The same 3-sigma upper convention is used for every trial; these are per-source
    # conditional bounds, not simultaneous confidence limits across the whole search.
    c, filters = continuum_upper(lam, 1.05 * photo_upper, filt)
    flo = 0.9 * max(0., line - 3*line_err)
    width = float(flo/(c*(1+z))) if c is not None and c > 0 else None
    return {"origin":origin, "aperture":aperture, "eligible":bool(eligible and width is not None),
            "r_upper_nanomaggy_before_cal_margin":float(photo_upper), "filters":filters,
            "line_flux":float(line), "line_flux_err":float(line_err),
            "line_flux_lower_cgs":float(flo), "continuum_upper_cgs_per_A":c,
            "ew0_lower_A":width}


def hsc_scale(hsc):
    units = str(hsc.get("unit","")).strip().lower().replace(" ","")
    if any(x in units for x in ("/sr","arcsec","steradian")):
        return None
    if "nanomagg" in units or "nmaggy" in units:
        return 1.
    if units in ("njy","njy/pixel","njy/pix"):
        return 1/3631.
    # Only interpret an AB zero point if pixel values are counts, or BUNIT is absent.
    if units not in ("","unreported","count","counts","adu","adu/pixel","adu/pix"):
        return None
    zero = num(hsc.get("fluxmag0",hsc.get("FLUXMAG0")))
    if zero and zero > 1e6:
        return 1e9/zero
    magzp = num(hsc.get("zeropoint",hsc.get("magzp",hsc.get("MAGZP"))))
    if magzp and 15 < magzp < 35:
        return 10**(0.4*(22.5-magzp))
    return None


def lan_test(rec, lan, legacy_upper=None):
    """Keep published optical and isophotal apertures explicit; reject r=30 sentinel."""
    if not lan:
        return None
    mag, merr = num(lan.get("hsc-r_mag")), num(lan.get("hsc-r_mag_err"))
    f, e = num(lan.get("flux_lya")), num(lan.get("flux_lya_err"))
    area = num(lan.get("area_r_iso_circ"))
    rad = math.sqrt(area/math.pi) if area and area>0 else None
    detail = {"catalog":"HLAN v0.3, Mentuch Cooper et al. 2026 ApJ 1000, 38, Table 2",
              "published_EW_A":num(lan.get("combined_eqw_rest_lya")), "r_AB":mag,"r_AB_err":merr,
              "isophotal_line_radius_arcsec":rad,"isophotal_line_flux_cgs":f,
              "isophotal_line_flux_err_cgs":e, "catalog_HSC_aperture":"2 arcsec in Table 2"}
    if not (f and e and 0 < e < f < 1e-9):
        detail["status"] = "missing or unverified physical line flux"
        return detail
    if legacy_upper is not None:
        detail["total_optical_test"] = bound(rec, "HLAN isophotal line / Legacy Tractor total r",
                                              LEGACY_R, legacy_upper, f, e,
                                              "isophotal line; total model optical flux")
    if mag is None or not (10 < mag < 29.99) or merr is None or not (0 < merr < .5):
        detail["status"] = "missing reliable HSC flux and uncertainty (r=30 is a sentinel)"
        return detail
    av = num(rec.get("av"),num(rec.get("catalog",{}).get("av")))
    if av is None or not 0 <= av < 5:
        detail["status"] = "missing foreground extinction"
        return detail
    fhi = 10**(0.4*(22.5-mag+3*merr)) * 10**(.4*av)
    matched = (rad is not None and rad <= 1.0)
    detail["catalog_aperture_test"] = bound(rec,"HLAN isophotal line / HSC aperture",
                                             HSC_R, fhi, f, e, "published apertures",
                                             eligible=matched)
    detail["status"] = "aperture-compatible" if matched else "different or unknown aperture; diagnostic only"
    return detail


def imaging_ew(rec, audit):
    if rec.get("channel") != "ew":
        return None
    line = num(rec.get("flux"))
    err = num(rec.get("flux_err"))
    tests = []
    legacy_upper = None
    legacy = audit.get("Legacy",{})
    matched = legacy.get("matches",[])
    if len(matched) == 1:
        obj = matched[0]
        val = num(obj.get("flux_r"))
        ivar = num(obj.get("flux_ivar_r"),0)
        trans = num(obj.get("mw_transmission_r"))
        if (num(obj.get("sep_arcsec"),99)<=1.25 and num(obj.get("maskbits"))==0
            and num(obj.get("nobs_r"),0)>0 and num(obj.get("fracmasked_r"),1)<.15
            and num(obj.get("fracflux_r"),99)<.5 and ivar>0 and trans and 0<trans<=1.001
            and val is not None):
            photo_err = 1/math.sqrt(ivar)
            legacy_upper = (max(0.,val)+3*photo_err)/trans
            test = bound(rec, "Legacy Tractor total r and PSF-extracted HETDEX line",
                         LEGACY_R, legacy_upper, line, err, "total model optical; PSF line")
            if test:
                test.update(flux_r_nanomaggy=val, error_r_nanomaggy=photo_err,
                            galactic_transmission_r=trans,ls_id=obj.get("ls_id"))
                tests.append(test)
    aperture_tests = []
    for key, filters in (("HSC", HSC_R), ("Legacy_forced", LEGACY_R)):
        hsc = audit.get(key,{})
        cube = audit.get("cube",{})
        center = cube.get("central_aperture",{})
        unit_scale = hsc_scale(hsc)
        radc, radh = num(cube.get("aperture_radius_arcsec")),num(hsc.get("aperture_radius_arcsec"))
        ra,dec = num(rec.get("ra")),num(rec.get("dec"))
        ra,dec = num(hsc.get("ra"),ra),num(hsc.get("dec"),dec)
        ra2,dec2 = num(cube.get("ra"),num(rec.get("ra_det"),ra)),num(cube.get("dec"),num(rec.get("dec_det"),dec))
        sep = None
        if None not in (ra,dec,ra2,dec2):
            sep = 3600*math.hypot(((ra2-ra+180)%360-180)*math.cos(math.radians(dec)),dec2-dec)
        f, e = num(hsc.get("flux")),num(hsc.get("flux_err"))
        fc, ec = num(center.get("value")),num(center.get("scatter"))
        if (unit_scale and f is not None and e and e>0 and fc and ec and ec>0 and
            num(hsc.get("blank_apertures"),0)>=16 and num(center.get("n"),0)>=16 and
            radc and radh and abs(radc-radh)<.05 and sep is not None and sep<.3):
            photo_upper = unit_scale*(max(0.,f)+3*e)
            test = bound(rec,key + " forced r / continuum-subtracted PDR1 cube",filters,
                         photo_upper,fc*1e-17,ec*1e-17,"matched 1.5 arcsec radius; observed fluxes")
            if test:
                test.update(flux_r_nanomaggy=f*unit_scale,error_r_nanomaggy=e*unit_scale,
                            center_offset_arcsec=sep, blank_apertures=hsc["blank_apertures"],
                            cube_control_windows=center["n"])
                aperture_tests.append(test)
    if aperture_tests:
        # Choose the more precise continuum image by blank-sky uncertainty before
        # considering target flux. A shallow nondetection does not constrain a deep one.
        preferred = min(aperture_tests, key=lambda t: t["error_r_nanomaggy"])
        for t in aperture_tests:
            if t is not preferred:
                t.update(eligible=False, comparison="shallower matched-aperture image")
            tests.append(t)
    lan = lan_test(rec,audit.get("LAN"),legacy_upper)
    if lan:
        tests.extend(x for k,x in lan.items() if k in ("total_optical_test","catalog_aperture_test") and x)
    eligible = [t for t in tests if t["eligible"]]
    return {"ew0_lower_A": min(t["ew0_lower_A"] for t in eligible) if eligible else None,
            "bounds":tests,"LAN_check":lan,"beta_interval":[-3,0],
            "rule":"Smallest eligible lower bound. Photometry +3 sigma and +5% scale, "
                   "line -3 sigma and -10% scale; band-integrated photon response. "
                   "At least 99% of the filter's AB weight must be redward of line+25 A. "
                   "Same-survey fluxes and LAN photometry can be correlated. Bounds are "
                   "conditional on the continuum-slope range and recorded apertures."}
