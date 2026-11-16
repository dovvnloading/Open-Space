"""Screen adjacent visits for paired changes in blue-green spectral curvature."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "cache"
BANDS = (413, 442, 465, 490, 510, 530, 560, 620, 665, 678)
# Bounds in degrees, inclusive; four nonoverlapping windows.
REGIONS = {
    "eastern_gulf": (24.1, 29.9, -88.95, -86.0),
    "west_florida_slope": (24.8, 29.9, -85.99, -82.5),
    "florida_straits": (24.1, 25.3, -82.49, -79.15),
    "florida_atlantic": (25.31, 29.9, -80.50, -79.15),
}


def read():
    keys = ["Rrs_" + str(b) for b in BANDS] + ["STRAYLIGHT_MASK", "time"]
    arrays = {k: [] for k in keys}
    paths = [CACHE / f"oci_{i:03}.nc" for i in range(0, 687, 36)]
    for p in paths:
        with netcdf_file(p, mmap=False) as ds:
            lat, lon = [ds.variables[k].data.copy() for k in ("latitude", "longitude")]
            for k in keys:
                arrays[k].append(ds.variables[k].data.copy())
    arrays = {k: np.concatenate(v) for k, v in arrays.items()}
    if len(arrays["time"]) != 172 or not np.all(np.diff(arrays["time"]) > 0):
        raise ValueError("Expected 172 ordered survey visits.")
    with netcdf_file(CACHE / "etopo.nc", mmap=False) as ds:
        y, x = [ds.variables[k].data for k in ("latitude", "longitude")]
        iy = np.abs(y[:, None] - lat).argmin(0)
        ix = np.abs(x[:, None] - lon).argmin(0)
        depth = -ds.variables["altitude"].data[np.ix_(iy, ix)]
    return arrays, lat, lon, depth


def window(lat, lon, bounds):
    s, n, w, e = bounds
    return (lat[:, None] >= s) & (lat[:, None] <= n) & (lon >= w) & (lon <= e)


def quality(a):
    waves = np.stack([a["Rrs_" + str(b)] for b in (442, 490, 560)])
    good = np.all(np.isfinite(waves) & (waves > 0) & (waves < 0.05), axis=0)
    green = a["Rrs_560"]
    red = a["Rrs_665"]
    return (
        (a["STRAYLIGHT_MASK"] == 0)
        & good
        & (green > 0.00005)
        & (red >= 0)
        & (red < 0.0015)
        & (red < 0.5 * green)
    )


def main():
    a, lat, lon, depth = read()
    dates = pd.to_datetime(a["time"], unit="s", utc=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        log_b, log_c, log_g = [np.log(a["Rrs_" + str(b)]) for b in (442, 490, 560)]
        curvature = log_c - (70 * log_b + 48 * log_g) / 118
    good = quality(a)
    visits = []
    pairs = []
    for reg, bounds in REGIONS.items():
        region = (depth >= 200) & window(lat, lon, bounds)
        total = int(region.sum())
        for t, d in enumerate(dates):
            use = good[t] & region
            row = dict(
                region=reg,
                date=d.date(),
                deep_points=total,
                usable_points=int(use.sum()),
            )
            if use.any():
                row["median_C"] = np.median(curvature[t, use])
                for wavelength in BANDS:
                    values = a[f"Rrs_{wavelength}"][t, use]
                    if np.isfinite(values).any():
                        row[f"median_Rrs_{wavelength}"] = np.nanmedian(values)
            visits.append(row)
            if not t:
                continue
            dt = (d - dates[t - 1]).days
            paired = region & good[t - 1] & good[t]
            n = int(paired.sum())
            row = dict(
                region=reg,
                start=dates[t - 1].date(),
                end=d.date(),
                days=dt,
                deep_points=total,
                paired_points=n,
            )
            row["coverage_pass"] = dt <= 8 and n >= 10 and n >= 0.15 * total
            if n:
                delta = (curvature[t] - curvature[t - 1])[paired]
                change = np.median(delta)
                row.update(
                    delta_C=change,
                    q25=np.quantile(delta, 0.25),
                    q75=np.quantile(delta, 0.75),
                    sign_agreement=np.mean(np.sign(delta) == np.sign(change)),
                    median_start_C=np.median(curvature[t - 1, paired]),
                    median_end_C=np.median(curvature[t, paired]),
                    median_lat=np.median(
                        np.broadcast_to(lat[:, None], region.shape)[paired]
                    ),
                    median_lon=np.median(
                        np.broadcast_to(lon[None, :], region.shape)[paired]
                    ),
                )
            pairs.append(row)
    pd.DataFrame(visits).to_csv(
        ROOT / "regional_visits.csv", index=False, float_format="%.8g"
    )
    out = pd.DataFrame(pairs)
    out.to_csv(ROOT / "regional_changes.csv", index=False, float_format="%.8g")
    passed = out[out.coverage_pass].copy()
    print(
        "Dates",
        len(dates),
        "grid",
        depth.shape,
        "passed",
        passed.groupby("region").size(),
    )
    print(
        passed.sort_values("delta_C", key=abs, ascending=False)
        .head(12)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
