"""Per-target Gaia transverse-speed evidence. Units: mas, mas/yr, km/s, pc."""
import csv
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
from scipy.optimize import brentq
from zero_point import zpt

ROOT = Path(__file__).resolve().parent
K = 4.74047
Z95 = NormalDist().inv_cdf(0.95)
# A sensitivity allowance in addition to formal errors (mas), not an estimate
# of an independent Gaia-wide random calibration term.
SYS = 0.020
KNOWN = {"5805243926609660032": "D6-1", "1798008584396457088": "D6-2", "2156908318076164224": "D6-3"}


def rows(name):
    with (ROOT / f"{name}.csv").open() as stream:
        return list(csv.DictReader(stream))


def getfloat(r, k):
    return float(r[k]) if r[k] else np.nan


def zp(r):
    return zpt.get_zpt([getfloat(r, "phot_g_mean_mag")], [getfloat(r, "nu_eff_used_in_astrometry")],
                       [getfloat(r, "pseudocolour")], [getfloat(r, "ecl_lat")], [int(r["astrometric_params_solved"])],
                       _warnings=False)[0]


def covariance(r):
    error = np.array([getfloat(r, c + "_error") for c in ("parallax", "pmra", "pmdec")])
    corrs = [getfloat(r, k) for k in ("parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr")]
    a, b, c = corrs
    corr = np.array([[1, a, b], [a, 1, c], [b, c, 1]])
    return corr * np.outer(error, error)


def quality(r):
    return (getfloat(r, "ruwe") < 1.4 and getfloat(r, "visibility_periods_used") >= 9
            and r["duplicated_source"] == "false" and getfloat(r, "ipd_frac_multi_peak") <= 2
            and getfloat(r, "ipd_gof_harmonic_amplitude") < 0.1)


def statistics(r, dist, fid, neighbours, ap):
    mu_vec = np.array([getfloat(r, k) for k in ("pmra", "pmdec")])
    mu = np.linalg.norm(mu_vec)
    pi = float(r["parallax"])
    zero_point = zp(r)
    assert np.isfinite(zero_point), r["source_id"]
    corrected = pi - zero_point
    cov = covariance(r)

    def score(v):
        # Positive iff 4.74047*mu/v - corrected parallax is positive.
        gradient = np.r_[-1., K * mu_vec / (v * mu)]
        sig = np.sqrt(gradient @ cov @ gradient + SYS**2)
        return (K * mu / v - corrected) / sig

    vlo = brentq(lambda v: score(v) - Z95, 10, 10000)
    dt = {k: float(dist[k]) * K * mu / 1000 for k in
          ("r_lo_geo", "r_med_geo", "r_hi_geo", "r_lo_photogeo", "r_med_photogeo", "r_hi_photogeo")}
    dr = math.pi / 180
    def sep(i):
        dec1, dec2 = float(r["dec"]) * dr, float(i["dec"]) * dr
        dlon = (float(r["ra"]) - float(i["ra"])) * dr
        h = math.sin((dec1-dec2)/2)**2 + math.cos(dec1)*math.cos(dec2)*math.sin(dlon/2)**2
        return math.degrees(2*math.asin(min(1, math.sqrt(h)))) * 3600
    near = [(sep(a), a) for a in neighbours if a["source_id"] != r["source_id"] and sep(a) <= 3]
    closest, companion = min(near, key=lambda x: x[0]) if near else (None, None)
    dmag = float(companion["phot_g_mean_mag"])-float(r["phot_g_mean_mag"]) if companion else None
    return dict(source_id=r["source_id"], name=KNOWN.get(r["source_id"], ""), raw_M_G=float(r["phot_g_mean_mag"])+5*math.log10(pi)-10,
                raw_vt=K*mu/pi, zp_mas=zero_point, corrected_vt=K*mu/corrected,
                z_margin_1000=score(1000), vt_lower95=vlo, **{"vt_"+k: v for k,v in dt.items()},
                quality_flags_pass=quality(r), fidelity_v2=fid["fidelity_v2"],
                neighbour_3arcsec=companion["source_id"] if companion else "", sep_arcsec=closest, delta_g=dmag,
                gspphot_teff=ap["teff_gspphot"], gspphot_logg=ap["logg_gspphot"],
                gspphot_d_pc=ap["distance_gspphot"],
                gspphot_model_vt=round(float(ap["distance_gspphot"])*K*mu/1000,2) if ap["distance_gspphot"] else "")


def main():
    zpt.load_tables()
    sources = rows("gaia")
    dist, fidelity, ap = ({r["source_id"]:r for r in rows(k)} for k in ("distance","fidelity","ap"))
    near = rows("gaia_3arcsec")
    results = [statistics(r, dist[r["source_id"]], fidelity[r["source_id"]], near, ap[r["source_id"]]) for r in sources]
    results.sort(key=lambda r:-r["vt_lower95"])
    with (ROOT / "kinematics.csv").open("w") as stream:
        w=csv.DictWriter(stream, fieldnames=list(results[0]), lineterminator="\n"); w.writeheader()
        for r in results:
            w.writerow({k: round(v, 6) if isinstance(v, (float, np.floating)) else v for k, v in r.items()})

if __name__ == "__main__":
    main()
