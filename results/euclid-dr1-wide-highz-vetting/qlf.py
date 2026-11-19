#!/usr/bin/env python3
"""Illustrative z>=7 quasar LF integration with explicit assumptions and citations.

Not a survey selection function, injection-recovery measurement, or forecast of
the number passing our strict dropout/NISP morphology cuts.
"""
import json
import math
from pathlib import Path

import numpy as np

OUT=Path(__file__).resolve().parent/"output"
H0=67.7
OM0=.31
C=299792.458
PHISTAR_GPC=1.35
MSTAR=-25.60
ALPHA=-1.20
BETA=-3.34
K=-.78
A_DEG2=1900.0
J_LIMIT_AB=23.
Z0,Z1=7.,9.


def expectation(area_deg2=A_DEG2,j_limit=J_LIMIT_AB,zmax=Z1):
    # Flat Lambda-CDM, homogeneous comoving volume.
    zs=np.linspace(0.,max(9.,zmax),18001)
    e=np.sqrt(OM0*(1+zs)**3+(1-OM0))
    dz=zs[1]-zs[0]
    dc=np.zeros_like(zs)
    dc[1:]=(C/H0)*np.cumsum((1/e[1:]+1/e[:-1])*dz*.5)
    zz=np.linspace(Z0,zmax,1001)
    e_interp=np.sqrt(OM0*(1+zz)**3+(1-OM0))
    r=np.interp(zz,zs,dc)
    dl=(1+zz)*r
    distance_modulus=5*np.log10(dl)+25
    # Flat AB f_nu approximation, with J chosen close to rest 1450 A.
    # No detailed filter-throughput, IGM, line or SED K-correction is claimed.
    ml= j_limit-distance_modulus+2.5*np.log10(1+zz)
    vol=(C/H0)*r**2/e_interp
    integrated_phi=[]
    for z,upper in zip(zz,ml):
        mag=np.linspace(-32.,float(upper),1601)
        dr=mag-MSTAR
        number=(PHISTAR_GPC*1e-9)*10**(K*(z-7)) / (
            10**(.4*(ALPHA+1)*dr)+10**(.4*(BETA+1)*dr))
        integrated_phi.append(np.trapz(number,mag))
    solid_angle=area_deg2*(math.pi/180)**2
    count=float(solid_angle*np.trapz(vol*np.asarray(integrated_phi),zz))
    return count, float(ml[0]),float(ml[-1])


def main():
    count,m0,m1=expectation()
    assert 17.5<count<20.5, "LF numerical integral outside independent benchmark"
    half=expectation(area_deg2=A_DEG2/2)[0]
    assert abs(half-count/2)<1e-6
    baseline=12*A_DEG2/3000.
    result={
      "method":"Double-power-law M1450 LF integrated over comoving volume and apparent J_AB limit",
      "lf_source":"Matsuoka et al. 2023 ApJL 949 L42 Table 2 standard fit; doi:10.3847/2041-8213/acd69f",
      "density_evolution":"k=-0.78 exponential-in-redshift normalization, held fixed when extrapolating beyond z=7",
      "k_source":"Matsuoka et al. 2023 Table 2 (fixed parameter), adopted also by Owl-z (2025)",
      "cosmology":{"H0_km_s_Mpc":H0,"Omega_m":OM0,"flat":True},
      "model":{"Mstar_AB":MSTAR,"alpha":ALPHA,"beta":BETA,
               "phi_star_Gpc_minus3_mag_minus1":PHISTAR_GPC,"k":K},
      "z_range":[Z0,Z1],"J_limit_AB":J_LIMIT_AB,
      "nominal_release_wide_area_deg2":A_DEG2,
      "predicted_J_le_23_z_7_to_9_nominal_area":round(count,3),
      "predicted_surface_density_deg_minus2":round(count/A_DEG2,6),
      "approx_M1450_at_limit_z7":round(m0,3),
      "approx_M1450_at_limit_z9":round(m1,3),
      "empirical_yang2026_z_gte_7_count":12,
      "empirical_yang2026_survey_area_approx_deg2":3000,
      "empirical_yang2026_scaled_nominal_area":round(baseline,3),
      "yang_reference":"Yang et al. 2026 A&A 711 A104, doi:10.1051/0004-6361/202658883",
      "not_a_completeness_corrected_strict_selection_expectation":True,
      "masked_release_area_measured":False,
      "limitations":[
        "Matsuoka 2023 LF fit at approximately z=7 extrapolated to z=9 with fixed k",
        "J treated as flat-fnu monochromatic AB flux near rest 1450 A without synthetic-filter SED and Ly-alpha transmission corrections",
        "Nominal 1900 deg2 not an exposure-weighted usable area",
        "J<=23 is a uniform threshold approximation; not the spatially varying detection or strict cut probability",
        "Different empirical and LF counts have selection, redshift, and flux-window differences; do not compare as discovery yield or a null result",
        "LF parameter covariance and cosmic variance not propagated"
      ]
    }
    OUT.mkdir(exist_ok=True)
    (OUT/"qlf_expectation.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("QLF_INTEGRAL",json.dumps({"count_nominal":round(count,3),"yang_scale":round(baseline,3)}),flush=True)


if __name__=="__main__":
    main()
