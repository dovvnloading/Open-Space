"""Fixed-sky aperture photometry and point-likeness on registered Stage 3 stamps."""

from pathlib import Path
import warnings, sys
import pandas as pd, numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates, gaussian_filter

warnings.simplefilter("ignore", FITSFixedWarning)
from catalog_screen import DEC
from collections import defaultdict

prefix = (
    "history_"
    if "--history" in sys.argv
    else (
        "blind_"
        if "--blind" in sys.argv
        else (
            "follow_"
            if "--follow" in sys.argv
            else ("onset_" if "--onsets" in sys.argv else "")
        )
    )
)
sites = pd.read_csv(
    "historical_sites.csv"
    if prefix == "history_"
    else (
        "blind_pairs.csv"
        if prefix == "blind_"
        else "onsets.csv" if prefix else "sites.csv"
    )
).set_index("site")
shifts = pd.read_csv("catalog_stats.csv").set_index(["obs", "filter"])
results = []
stamps = defaultdict(dict)
fluxes = []
# The 0.03 arcsec sky grid represents surface brightness. Apertures are NOT point-source corrected.
scale = 0.03
gy, gx = np.mgrid[-44:45, -44:45] * scale
rr = np.hypot(gx, gy)
pix_njy = (np.deg2rad(scale / 3600)) ** 2 * 1e15
ap = rr <= 0.15
ann = (rr >= 0.27) & (rr <= 0.42)
offsets = [
    (x, y)
    for x in np.arange(-0.9, 0.901, 0.3)
    for y in np.arange(-0.9, 0.901, 0.3)
    if 0.45 <= np.hypot(x, y) <= 1.06
]
masks = [
    (
        np.hypot(gx - x, gy - y) <= 0.15,
        (np.hypot(gx - x, gy - y) >= 0.27) & (np.hypot(gx - x, gy - y) <= 0.42),
    )
    for x, y in offsets
]


def flux(image, aper, ring):
    if not np.isfinite(image[aper]).all() or np.isfinite(image[ring]).mean() < 0.9:
        return np.nan
    return np.sum(image[aper] - np.nanmedian(image[ring])) * pix_njy


def measure_scatter(values):
    good = np.isfinite(values)
    return (
        1.4826 * np.nanmedian(np.abs(values - np.nanmedian(values)))
        if sum(good) >= 15
        else np.nan
    )


for path in sorted(Path(prefix + "stamps").glob("*.npz")):
    pid, obs, band = map(int, path.stem.split("-"))
    fil = f"F{band:03d}W"
    p = np.load(path)
    h = fits.Header.fromstring(str(p["header"]), sep="\n")
    w = WCS(h, naxis=2)
    date = str(p["date"])
    for i, site in enumerate(p["site"]):
        # Correct the 2025 October WCS to the new-epoch reference using common objects.
        s = shifts.loc[obs, fil] if pid == 8060 and obs <= 2 else None
        dx = s.offset_east_arcsec if s is not None else 0
        dy = s.offset_north_arcsec if s is not None else 0
        row = sites.loc[site]
        cos = np.cos(np.deg2rad(row.dec))
        x, y = w.all_world2pix(
            row.ra + (gx - dx) / (3600 * cos), row.dec + (gy - dy) / 3600, 0
        )
        im = map_coordinates(
            p["data"][i],
            [y - p["y0"][i], x - p["x0"][i]],
            order=1,
            mode="constant",
            cval=np.nan,
        )
        fb = flux(im, ap, ann)
        blanks = np.array([flux(im, a, b) for a, b in masks])
        e = measure_scatter(blanks)
        if not np.isfinite(fb) or not np.isfinite(e) or e <= 0:
            continue
        fluxes.append(
            dict(
                site=site,
                pid=pid,
                obs=obs,
                filter=fil,
                utc=date,
                aperture_njy=fb,
                local_scatter_njy=e,
            )
        )
        if obs in (1, 2) and pid in [8060, 12577]:
            stamps[site, obs, fil][pid] = im
for (site, obs, fil), pair in stamps.items():
    if len(pair) != 2:
        continue
    diff = pair[12577] - pair[8060]
    df = flux(diff, ap, ann)
    blanks = np.array([flux(diff, a, b) for a, b in masks])
    err = measure_scatter(blanks)
    if not np.isfinite(df) or not np.isfinite(err) or err <= 0:
        continue
    sign = np.sign(df)
    z = sign * (diff - np.nanmedian(diff[ann]))
    yy, xx = np.nonzero(rr <= 0.12)
    best = np.nanargmax(gaussian_filter(np.nan_to_num(z), 0.85)[yy, xx])
    x0 = gx[yy[best], xx[best]]
    y0 = gy[yy[best], xx[best]]
    r = np.hypot(gx - x0, gy - y0)
    w = r < 0.20
    sig = (int(fil[1:4]) / 100 * 0.034 + 0.025) / 2.355
    template = np.exp(-0.5 * (r / sig) ** 2)
    corr = np.corrcoef(template[w], z[w])[0, 1] if np.isfinite(z[w]).all() else np.nan
    near = rr < 0.27
    pos = np.sum(np.maximum(z[near], 0))
    neg = np.sum(np.maximum(-z[near], 0))
    results.append(
        dict(
            site=site,
            obs=obs,
            filter=fil,
            new_minus_old_ap_njy=df,
            local_scatter_njy=err,
            significance=df / err,
            peak_offset_arcsec=np.hypot(x0, y0),
            gaussian_shape_corr=corr,
            negative_positive_ratio=neg / pos if pos else np.nan,
            n_control_apertures=np.isfinite(blanks).sum(),
        )
    )
pd.DataFrame(fluxes).sort_values(["site", "utc", "filter"]).to_csv(
    prefix + "forced_lightcurves.csv", index=False
)
if prefix == "history_":
    sys.exit(0)
r = pd.DataFrame(results).sort_values(["site", "obs", "filter"])
r.to_csv(prefix + "image_checks.csv", index=False)
g = r[
    (abs(r.significance) >= 5)
    & (r.peak_offset_arcsec <= 0.10)
    & (r.gaussian_shape_corr >= 0.65)
    & (r.negative_positive_ratio < 0.65)
]
g = g.assign(sign=np.sign(g.significance))
passed = g.groupby(["site", "sign"])["filter"].nunique()
sel = passed[passed >= 2]
print(len(r), "measurements, multi-band morphology survivors:", sel.to_dict())
print(
    g[g.site.isin(sel.index.get_level_values(0))][
        ["site", "filter", "significance", "gaussian_shape_corr"]
    ]
    .round(2)
    .to_string(index=False)
)
