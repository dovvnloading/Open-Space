"""Toy two-image singular-isothermal-sphere (SIS) check, NOT a lens measurement.

Lemon et al. measured source z=1.62, but not a secure lens galaxy redshift.
Sweep *assumed* lens z. This toy deliberately ignores the three observed lens
galaxies, possible blended extra image, extinction, shear and microlensing.
"""
from math import pi, sqrt
from scipy.constants import c, G, parsec
from scipy.integrate import quad

SOURCE_Z = 1.618
SEP_ARCSEC = 5.83
GMAG_A, GMAG_B = 17.407343, 19.373940
THETA_E = (SEP_ARCSEC / 2) * pi / (180 * 3600)
MPC = parsec * 1e6
H0_SI = 70 * 1000 / MPC
OM = .3
ZS = SOURCE_Z


def comoving(z):
    return c / H0_SI * quad(lambda x: 1 / sqrt(OM*(1+x)**3 + 1-OM), 0, z)[0]


def sis(zlens, delta_days=210):
    chi_s = comoving(ZS)
    chi_l = comoving(zlens)
    dl, ds, dls = chi_l/(1+zlens), chi_s/(1+ZS), (chi_s-chi_l)/(1+ZS)
    d_delay = (1+zlens)*dl*ds/dls
    max_delay_days = (d_delay/c) * 2*THETA_E**2 / 86400
    beta_by_theta_e = delta_days/max_delay_days
    sigma_v = c * sqrt((THETA_E/(4*pi))*ds/dls)/1000
    mass = c**2/(4*G) * dl*ds/dls*THETA_E**2/1.98847e30
    ratio = (1 + beta_by_theta_e)/(1 - beta_by_theta_e) if beta_by_theta_e < 1 else float('nan')
    return (max_delay_days, beta_by_theta_e, sigma_v, mass, ratio)


if __name__ == "__main__":
    obs_ratio = 10**(.4*(GMAG_B-GMAG_A))
    print("published_source_z", SOURCE_Z, "lens_redshift_unmeasured", True)
    print("image_sep_arcsec", SEP_ARCSEC, "observed_G_flux_ratio", round(obs_ratio, 3))
    print("ASSUMED_z_l max_SIS_delay_days beta_over_thetaE SIS_sigma_kms lens_Msun SIS_flux_ratio")
    for zl in (.2,.3,.4,.5,.6,.8,1.0):
        print(zl, *[round(v,3) for v in sis(zl)])
