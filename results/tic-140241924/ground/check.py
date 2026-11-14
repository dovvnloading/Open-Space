"""Test timing curvature in ASAS-SN with each binary's linear ephemeris free."""
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
CLOCKS = [("A primary", 1.8035, 1544.5094, .014),
          ("A secondary", 1.8035, 1544.5094 + 1.8035*.5042, .015),
          ("B primary", 4.3829, 1561.2557, .007)]


def prepare():
    ground = pd.read_csv(ROOT / "photometry.csv.gz")
    ground = ground[(ground.quality == "G") & ground.mag_err.between(.001, .2)].sort_values("jd").copy()
    jd = Time(ground.jd.to_numpy(), format="jd", scale="utc")
    ra, dec = np.deg2rad([126.771492, -50.765804])
    direction = np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])
    earth = get_body_barycentric("earth", jd.tdb).xyz.value
    ground["btjd"] = jd.tdb.jd - 2457000 + direction @ earth * 499.00478 / 86400
    seasons = pd.read_csv(ROOT.parent / "seasons.csv").query('aperture=="box" and window==0.17')
    first = seasons[seasons.clock == "A primary"].mean_btjd
    ground = ground[ground.btjd.between(first.min(), first.max())].copy()
    time = ground.btjd.to_numpy()
    for clock, _, _, _ in CLOCKS:
        g = seasons[seasons.clock == clock]
        t, y = g.mean_btjd.to_numpy(), g.delta.to_numpy()
        endpoint_line = y[0] + (y[-1]-y[0])*(time-t[0])/(t[-1]-t[0])
        ground[clock + " line"] = endpoint_line
        ground[clock + " curve"] = np.interp(time, t, y) - endpoint_line
    return ground


def fit(ground, scales):
    t, mag = ground.btjd.to_numpy(), ground.mag.to_numpy()
    error = np.hypot(.02, ground.mag_err.to_numpy())
    baseline_groups = ground.camera.astype(str) + "_" + (t//365).astype(int).astype(str)
    mat = pd.get_dummies(baseline_groups, dtype=float).to_numpy()
    n = mat.shape[1]
    time = (t-t.min())/(t.max()-t.min()) - .5

    def residual(params):
        model = mat @ params[:n]
        shifts = params[n+3:n+7]
        for j, (clock, period, epoch, width) in enumerate(CLOCKS):
            binary = int(j == 2)
            linear = shifts[2*binary] + shifts[2*binary+1]*time
            lag = ground[clock + " line"].to_numpy() + scales[binary]*ground[clock + " curve"].to_numpy() + linear
            phase = ((t-epoch-lag/1440)/period + .5) % 1 - .5
            model += params[n+j] * np.exp(-.5*(phase/width)**2)
        return (model-mag)/error

    init = np.r_[np.linalg.lstsq(mat, mag, rcond=None)[0], [.2,.08,.1], np.zeros(4)]
    lower = np.r_[[13.]*n, [0.]*3, [-180.]*4]
    upper = np.r_[[18.]*n, [1.5]*3, [180.]*4]
    result = least_squares(residual, init, bounds=(lower,upper), loss="huber", f_scale=2.5)
    if not result.success:
        raise RuntimeError(result.message)
    return dict(loss=result.cost, A_primary_mag=result.x[n], A_secondary_mag=result.x[n+1],
                B_primary_mag=result.x[n+2])


if __name__ == "__main__":
    data = prepare()
    samples = {"all": data}
    samples.update({f"omit camera {c}": data[data.camera != c] for c in sorted(data.camera.unique())})
    year = (data.btjd//365).astype(int)
    samples.update({f"omit 365-day bin {i}": data[year != i] for i in sorted(year.unique()) if (year == i).sum() > 40})
    rows = []
    for group, sample in samples.items():
        for a, b in [(0,0),(1,1),(1,0),(0,1),(-1,-1)]:
            rows.append(dict(subset=group, observations=len(sample), A_scale=a, B_scale=b,
                             **fit(sample,[a,b])))
    pd.DataFrame(rows).to_csv(ROOT / "fits.csv", index=False, float_format="%.6f")
