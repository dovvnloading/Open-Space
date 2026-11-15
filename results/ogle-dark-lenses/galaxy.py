"""Galactic prior for lens mass and distance from (tE, piE) samples.

For each posterior sample the unmeasured relative proper motion and source
distance are drawn and weighted with the event-rate expression of Batista et
al. (2011, eq. 18):

    weight = nu(D_L) f(mu) [g(M) M] D_L^4 mu^4 tE / piE * nu_s(D_S)

with theta_E = mu tE, M = theta_E / (kappa piE), pi_rel = theta_E piE.
The lens density uses a double-exponential disc and the Han & Gould (2003)
bar; kinematics are Gaussians about disc rotation or a non-rotating bulge.
Quantities are heliocentric: the geocentric fit is corrected with Earth's
projected velocity at the reference time.
"""
import numpy as np

KAPPA = 8.144            # mas / Msun
R0 = 8.2                 # kpc
V0 = 233.0               # km/s, circular speed at the Sun
SUN = (11.1, 12.24, 7.25)  # km/s, solar peculiar motion (U, V, W)
KMS_PER_MASYR_KPC = 4.74047

# Main-sequence absolute I magnitude against mass (Pecaut & Mamajek 2013 dwarf table).
MS_MASS = np.array([0.10, 0.16, 0.23, 0.37, 0.50, 0.57, 0.70, 0.88, 1.00, 1.33, 2.3, 4.7, 7.3, 17.7, 40.0])
MS_MI = np.array([12.5, 11.0, 10.0, 8.9, 7.9, 7.15, 6.03, 5.05, 4.13, 2.9, 1.1, -0.75, -1.45, -3.2, -4.6])


def galactic(ra, dec):
    """Galactic longitude and latitude (rad) plus the position angles (rad,
    north through east) of increasing longitude and of the north Galactic pole."""
    ra_gp, dec_gp, l_ncp = np.radians(192.85948), np.radians(27.12825), np.radians(122.93192)
    b = np.arcsin(np.sin(dec_gp) * np.sin(dec) + np.cos(dec_gp) * np.cos(dec) * np.cos(ra - ra_gp))
    y = np.cos(dec) * np.sin(ra - ra_gp)
    x = np.cos(dec_gp) * np.sin(dec) - np.sin(dec_gp) * np.cos(dec) * np.cos(ra - ra_gp)
    l = (l_ncp - np.arctan2(y, x)) % (2 * np.pi)
    angle_ngp = np.arctan2(np.cos(dec_gp) * np.sin(ra_gp - ra),
                           np.sin(dec_gp) * np.cos(dec) - np.cos(dec_gp) * np.sin(dec) * np.cos(ra_gp - ra))
    return l, b, angle_ngp + np.pi / 2, angle_ngp


def galactocentric(dist, l, b):
    """Cartesian offsets from the Galactic centre (kpc): x along the Sun-centre
    line (positive beyond the centre), y towards l = 90 deg, z to the north pole."""
    return dist * np.cos(b) * np.cos(l) - R0, dist * np.cos(b) * np.sin(l), dist * np.sin(b)


def density(dist, l, b):
    """Disc and bar mass densities (Msun/pc^3) at heliocentric distance (kpc)."""
    x, y, z = galactocentric(dist, l, b)
    radius = np.hypot(x, y)
    disc = 0.05 * np.exp(-(radius - R0) / 2.6 - np.abs(z) / 0.3)
    angle = np.radians(25.0)        # near end of the bar at positive longitude
    xb = -x * np.cos(angle) + y * np.sin(angle)
    yb = x * np.sin(angle) + y * np.cos(angle)
    rs = (((xb / 1.58) ** 2 + (yb / 0.62) ** 2) ** 2 + (z / 0.43) ** 4) ** 0.25
    bar = 1.23 * np.exp(-0.5 * rs ** 2)
    return disc, bar


def mean_velocity(dist, l, b, rotating):
    """Mean velocity of a population relative to the Sun along increasing l and
    b (km/s). `rotating` scales a flat rotation curve (1 disc, 0 bulge)."""
    x, y, _ = galactocentric(dist, l, b)
    radius = np.maximum(np.hypot(x, y), 1e-3)
    vx = V0 * rotating * y / radius - SUN[0]
    vy = -V0 * rotating * x / radius - (V0 + SUN[1])
    vz = -SUN[2]
    return (-vx * np.sin(l) + vy * np.cos(l),
            -vx * np.sin(b) * np.cos(l) - vy * np.sin(b) * np.sin(l) + vz * np.cos(b))


def ms_apparent_i(mass, dist, a_i_total, b, d_source):
    """I magnitude of a main-sequence lens behind a share of the total extinction."""
    scale = 0.1 / max(abs(np.sin(b)), 1e-3)      # kpc along the sight line per dust scale height
    share = (1 - np.exp(-dist / scale)) / (1 - np.exp(-d_source / scale))
    return np.interp(np.log10(mass), np.log10(MS_MASS), MS_MI) + 5 * np.log10(dist * 100) + a_i_total * share


def lens_posterior(te, pi_n, pi_e, blend_flux, ra, dec, earth_velocity_ne, a_i_total=1.5,
                   draws=40, mass_slope=-1.0, seed=1, blend_zero_point=18.0):
    """Weighted samples of lens mass and distance.

    te (days), pi_n, pi_e: geocentric posterior samples.
    earth_velocity_ne: Earth's velocity projected north, east (au/yr) at t0par.
    mass_slope: g(M) ~ M^slope. Returns a dict of arrays with weights.
    """
    rng = np.random.default_rng(seed)
    l, b, angle_l, angle_b = galactic(ra, dec)
    n = len(te)
    te, pi_n, pi_e, blend_flux = (np.repeat(np.asarray(a, float), draws) for a in (te, pi_n, pi_e, blend_flux))
    pie = np.hypot(pi_n, pi_e)
    mu_geo = rng.uniform(0.0, 25.0, n * draws)                # mas/yr
    d_s = rng.uniform(4.0, 12.0, n * draws)
    theta = mu_geo * te / 365.25                              # mas
    mass = theta / (KAPPA * pie)
    pi_rel = theta * pie
    d_l = 1.0 / (pi_rel + 1.0 / d_s)
    # heliocentric relative proper motion vector (north, east)
    mu_n = mu_geo * pi_n / pie + earth_velocity_ne[0] * pi_rel
    mu_e = mu_geo * pi_e / pie + earth_velocity_ne[1] * pi_rel
    mu_l = mu_n * np.cos(angle_l) + mu_e * np.sin(angle_l)
    mu_b = mu_n * np.cos(angle_b) + mu_e * np.sin(angle_b)
    disc_l, bar_l = density(d_l, l, b)
    disc_s, bar_s = density(d_s, l, b)
    source_mean = mean_velocity(d_s, l, b, 0.0)
    source_sigma = 100.0
    weight = np.zeros(n * draws)
    part = {}
    for name, rho, rotating, sig_l, sig_b in (("disc", disc_l, 1.0, 34.0, 18.0), ("bulge", bar_l, 0.0, 100.0, 100.0)):
        lens_mean = mean_velocity(d_l, l, b, rotating)
        mean_l = (lens_mean[0] / d_l - source_mean[0] / d_s) / KMS_PER_MASYR_KPC
        mean_b = (lens_mean[1] / d_l - source_mean[1] / d_s) / KMS_PER_MASYR_KPC
        var_l = (sig_l / d_l) ** 2 / KMS_PER_MASYR_KPC ** 2 + (source_sigma / d_s) ** 2 / KMS_PER_MASYR_KPC ** 2
        var_b = (sig_b / d_l) ** 2 / KMS_PER_MASYR_KPC ** 2 + (source_sigma / d_s) ** 2 / KMS_PER_MASYR_KPC ** 2
        f_mu = np.exp(-0.5 * ((mu_l - mean_l) ** 2 / var_l + (mu_b - mean_b) ** 2 / var_b)) / (
            2 * np.pi * np.sqrt(var_l * var_b))
        part[name] = rho * f_mu
        weight += part[name]
    mu_hel = np.hypot(mu_n, mu_e)
    weight *= mass ** (mass_slope + 1) * d_l ** 4 * mu_hel ** 4 * te / pie * (disc_s + bar_s)
    weight[~np.isfinite(weight) | (d_l <= 0.05) | (d_l >= d_s)] = 0.0
    lens_i = ms_apparent_i(np.clip(mass, 0.08, 60), np.clip(d_l, 0.05, None), a_i_total, b, d_s)
    blend_i = np.where(blend_flux > 0, blend_zero_point - 2.5 * np.log10(np.maximum(blend_flux, 1e-9)), 99.0)
    return {"mass": mass, "d_l": d_l, "d_s": d_s, "mu_hel": mu_hel, "theta_e": theta, "weight": weight,
            "ms_lens_i": lens_i, "blend_i": blend_i, "disc_fraction": part["disc"] / (part["disc"] + part["bulge"] + 1e-300)}


def quantiles(values, weights, q=(0.16, 0.5, 0.84)):
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    cumulative /= cumulative[-1]
    return np.interp(q, cumulative, values[order])
