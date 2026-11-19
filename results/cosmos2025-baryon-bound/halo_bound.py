#!/usr/bin/env python3
"""Survey-volume extreme-halo mass expectations under two HMF prescriptions.

Compute the distribution of the most massive halo in each *specified bin*:
    lambda(>M) = integral_z dV_com/dz integral_lnM^infty (dn/dlnM) dlnM
    P(M_max <= M) = exp[-lambda(>M)]   (Poisson, neglecting sample variance).
Upper bounds on stellar mass under maximally efficient baryon conversion:
    M_star <= f_b M_halo,  f_b = 0.16.
This is NOT an absolute hard upper mass, nor the likelihood of a given
galaxy, nor an Eddington- or cosmic-variance-corrected survey-level test.

Explicit mass definitions:
  ST99 from Colossus: FOF (Friends of Friends), not 200m
  Tinker08 from Colossus: M_200mean (200 times mean matter density)
Incompatibility of the definitions is a systematic, not a statistical error.
Tinker08 high-z application is an extrapolation beyond calibration.
"""
import argparse
import datetime
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq
from scipy.integrate import quad
from colossus.cosmology import cosmology
from colossus.lss import mass_function

FB = 0.16
AREA_DEG2 = 0.54  # approximate NIRCam coverage; not a depth/completeness mask
BINS = ((6., 7.), (7., 8.), (8., 9.), (9., 10.), (10., 12.), (12., 15.))
H0 = 67.32
OM = 0.3158
OB = 0.0493
S8 = 0.812
NS = 0.96605
C = 299792.458
MASS_GRID = np.logspace(8.0, 17.0, 500)  # Msun physical
Z_GRID_POINTS = 31  # full volume-integrated HMF, not z-midpoint only


def hz(z):
    # Flat LCDM; radiation adds at most ~0.3% to H at z=15. Document this approximation.
    return H0 * np.sqrt(OM * (1 + z)**3 + (1 - OM))


def dc(z):
    return C * quad(lambda zz: 1. / hz(zz), 0., z, epsrel=1e-8)[0]


def dVdz(z):
    area_sr = AREA_DEG2 * (np.pi / 180)**2
    return area_sr * C / hz(z) * dc(z)**2  # Mpc^3 per unit z


def threshold(masses, lam, quantile):
    # lambda(m) monotonically decreases with increasing M; CDF(max <= M) = exp(-lambda).
    required = -math.log(quantile)
    if lam[0] < required or lam[-1] > required:
        return None
    return float(10**np.interp(np.log(required), np.log(lam[::-1].clip(1e-250)), np.log(masses[::-1])))


def calculate():
    cosmo = cosmology.setCosmology("open_space_bk2023", params={
        "flat": True, "H0": H0, "Om0": OM, "Ob0": OB, "sigma8": S8, "ns": NS
    })
    h = H0 / 100
    models = [("sheth99", "fof"), ("tinker08", "200m")]
    result = {
        "computed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "area_deg2_nominal": AREA_DEG2, "f_baryon": FB,
        "cosmology": {"H0": H0, "Om0": OM, "Ob0": OB, "sigma8": S8, "ns": NS,
                      "flat": True, "linear_transfer": "Colossus default Eisenstein & Hu 1998"},
        "max_distribution": "P(M_max <= M) = exp(-V integral_{M}^{infty} dn/dlnM dlnM); Poisson approximation",
        "model_uncertainties": ["MFOF versus M200m are not the same halo mass definition",
                                "Tinker08 beyond calibration redshifts",
                                "no cosmic variance, baryon depletion, completeness, or lensing corrections",
                                "all nominal NIRCam area assumed valid in every bin; real selection varies"],
        "bins": []
    }
    for z1,z2 in BINS:
        zz = np.linspace(z1, z2, Z_GRID_POINTS)
        vol = float(trapezoid([dVdz(z) for z in zz], zz))
        b = {"z_min_exclusive": z1, "z_max_inclusive": z2, "volume_mpc3": vol,
             "models": {}}
        for model, mdef in models:
            m_h = MASS_GRID * h
            cumul = []
            for z in zz:
                # Colossus: dn/d ln(M) in comoving (Mpc/h)^-3, M in Msun/h.
                dn = mass_function.massFunction(m_h, z, mdef=mdef, model=model, q_out="dndlnM")
                if not np.all(np.isfinite(dn)) or np.any(dn < 0):
                    raise ValueError(f"bad halo number density: {model} at z={z}")
                # descending integration produces negative integral; flip sign and reverse
                n_gt = -cumulative_trapezoid(dn[::-1], np.log(m_h[::-1]), initial=0.)[::-1] * h**3
                cumul.append(n_gt * dVdz(z))
            # Expected number lambda across the complete redshift shell for each M threshold.
            lam = trapezoid(np.asarray(cumul), zz, axis=0)
            if np.any(np.diff(lam) > 1e-5) or np.any(lam < -1e-8):
                raise ValueError("non-monotone integrated abundance")
            quant = {f"p{int(q*1000):03d}": threshold(MASS_GRID, lam, q)
                     for q in (0.5, 0.95, 0.99, 0.999)}
            m1 = float(10**np.interp(np.log(1.), np.log(lam[::-1].clip(1e-250)), np.log(MASS_GRID[::-1])))
            b["models"][model] = {
                "mass_definition": mdef,
                "M_halo_msun_expected_one": m1,
                "log10_M_halo_expected_one": math.log10(m1),
                "log10_M_star_baryon_expected_one": math.log10(FB * m1),
                "max_halo_quantiles_msun": quant,
                "max_stellar_baryon_quantiles_log10": {k: math.log10(FB*x) if x is not None else None
                                                       for k,x in quant.items()},
                "lambda_at_halo_mass_1e12": float(np.interp(12., np.log10(MASS_GRID), lam)),
                "lambda_at_halo_mass_1e13": float(np.interp(13., np.log10(MASS_GRID), lam)),
            }
        result["bins"].append(b)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args=p.parse_args()
    result=calculate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"output": str(args.out), "bins": [
        {"range":(v["z_min_exclusive"],v["z_max_inclusive"]),
         "volume_mpc3": v["volume_mpc3"],
         "ST_logMstar_p95": v["models"]["sheth99"]["max_stellar_baryon_quantiles_log10"]["p950"]}
        for v in result["bins"]]},indent=2))


if __name__ == "__main__":
    main()
