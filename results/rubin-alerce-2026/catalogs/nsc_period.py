"""Independent NSC DR2 check at the fixed Rubin candidate period; no refit of frequency."""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

HERE = Path(__file__).resolve().parent
iers.conf.auto_download = False


def test():
    d = pd.read_csv(HERE / "nsc-epochs.csv", dtype={"rubin_oid": str})
    d = d.loc[
        (d.rubin_oid == "313761042359910519")
        & (d["flags"] == 0)
        & (d.magerr_auto > 0)
        & (d.magerr_auto < 1)
    ].copy()
    # The NSC exposure MJDs are used as archived. Unknown exposure start/mean
    # convention and band-dependent exposure lengths limit the phase test.
    t = Time(
        d.mjd.values,
        format="mjd",
        scale="utc",
        location=EarthLocation.from_geodetic(-70.8065, -30.169, 2200 * u.m),
    )
    target = SkyCoord(61.057494, -49.129764, unit="deg")
    bjd = (t.tdb + t.light_travel_time(target)).jd
    sigma = d.magerr_auto.values
    mag = d.mag_auto.values
    offsets = np.column_stack(
        [d["filter"].values == b for b in sorted(d["filter"].unique())]
    )

    def residual(A):
        x, *_ = np.linalg.lstsq(A / sigma[:, None], mag / sigma, rcond=None)
        return np.sum(((mag - A @ x) / sigma) ** 2), x

    chi0, _ = residual(offsets)
    period = 5.343047965998457
    theta = (bjd - bjd.min()) * (2 * np.pi * 1440 / period)
    chi1, coefs = residual(np.c_[offsets, np.sin(theta), np.cos(theta)])
    return dict(
        n=len(d),
        fixed_period_minutes=period,
        chi_square_constant=chi0,
        chi_square_sinusoid=chi1,
        delta_chi_square=chi0 - chi1,
        fitted_semia_magnitudes=np.hypot(*coefs[-2:]),
    )


if __name__ == "__main__":
    print(json.dumps(test(), indent=2))
