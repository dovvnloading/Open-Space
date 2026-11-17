"""Conditional injection recovery using Gaia's random_index as a sampling frame."""
import csv
from collections import defaultdict

import numpy as np
from zero_point import zpt

from analysis import K, ROOT, SYS, Z95, covariance, getfloat, quality, rows, zp

SEED = 260611293
NTRIAL = 256
DISTANCES_KPC = (1, 2, 4, 8)
SPEEDS = (500, 1000, 1500, 2500)
BINS = ((16, 17), (17, 18), (18, 19), (19, 20.5))


def simulate(row, rng):
    gmag = getfloat(row, "phot_g_mean_mag")
    mag_bin = next((lo, hi) for lo, hi in BINS if lo <= gmag < hi)
    valid = int(row["astrometric_params_solved"]) in (31, 95)
    zp_mas = zp(row) if valid else np.nan
    calibrated = np.isfinite(zp_mas)
    if valid and not calibrated:
        # These raw triggers still exist, but lack an in-domain individual ZP.
        zp_mas = -0.017
    if valid:
        cov = covariance(row)
        # Error ellipsoid and orientation are resampled once per template;
        # paired random numbers reduce noise in comparisons between conditions.
        noise = rng.multivariate_normal(np.zeros(3), cov, NTRIAL, check_valid="raise")
        azimuth = rng.uniform(0, 2*np.pi, NTRIAL)
        tangential = np.sqrt(1 - rng.uniform(-1, 1, NTRIAL)**2)
        unit_ra, unit_dec = np.cos(azimuth)*tangential, np.sin(azimuth)*tangential
        q_flag = quality(row)
    for distance in DISTANCES_KPC:
        for speed in SPEEDS:
            key = (*mag_bin, distance, speed)
            if not valid:
                # A two-parameter astrometric solution cannot pass the raw trigger.
                yield key, (0.,0.)
                continue
            parallax = 1/distance + zp_mas + noise[:, 0]
            mua = speed/(K*distance)*unit_ra + noise[:, 1]
            mud = speed/(K*distance)*unit_dec + noise[:, 2]
            mu = np.hypot(mua, mud)
            # Rearrangement of M_G > 3.5 avoids log of negative observations.
            pmin_cmd = 10**((3.5+10-gmag)/5)
            raw = ((mu > 35) & (parallax > 3*np.sqrt(cov[0,0])) &
                   (K*mu > 1000*parallax) & (parallax > pmin_cmd))
            da, dd = K/1000*mua/mu, K/1000*mud/mu
            variance = (cov[0,0] + da**2*cov[1,1] + dd**2*cov[2,2] - 2*da*cov[0,1]
                        - 2*dd*cov[0,2] + 2*da*dd*cov[1,2] + SYS**2)
            score = (K*mu/1000 - (parallax-zp_mas))/np.sqrt(variance)
            strong = raw & calibrated & q_flag & (score > Z95)
            yield key, (raw.mean(),strong.mean())


def main():
    rng = np.random.default_rng(SEED)
    zpt.load_tables()
    rates = defaultdict(list)
    for row in sorted(rows("templates"), key=lambda r:int(r["random_index"])):
        for key, probabilities in simulate(row, rng):
            rates[key].append(probabilities)
    with (ROOT / "injection.csv").open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["g_lo","g_hi","distance_kpc","true_3d_kms","templates","trials",
                         "raw_fraction","strong_fraction","raw_template_se","strong_template_se"])
        for key, samples in sorted(rates.items()):
            a = np.array(samples)
            # Template-level independent units, not one independent real star per draw.
            means = a.mean(axis=0)
            errors = a.std(axis=0, ddof=1)/np.sqrt(len(a))
            writer.writerow([*key,len(a),len(a)*NTRIAL,*np.round(means, 6),*np.round(errors, 6)])

if __name__ == "__main__":
    main()
