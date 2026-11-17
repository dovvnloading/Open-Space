"""Summarize provisional nuclear variations with matched rest-band aperture difference magnitudes."""

import numpy as np, pandas as pd
from astropy.cosmology import Planck18
from pathlib import Path
from scipy.spatial import cKDTree

sites = pd.read_csv("sites.csv").set_index("site")
images = pd.read_csv("image_checks.csv")
hosts = pd.read_csv("host_neighbors.csv").query("neighbor_rank==1").set_index("site")
c = images[
    (abs(images.significance) >= 5)
    & (images.gaussian_shape_corr >= 0.65)
    & (images.peak_offset_arcsec <= 0.10)
    & (images.negative_positive_ratio < 0.65)
].copy()
c["sign"] = np.sign(c.significance)
n = c.groupby(["site", "sign"])["filter"].nunique()
target = n[n >= 2].index.get_level_values("site")
flux = pd.read_csv("dither_fluxes.csv")
records = []
fn = lambda m, z: m - Planck18.distmod(z).value + 2.5 * np.log10(1 + z)
for site in target:
    h = hosts.loc[site]
    i = images[(images.site == site) & (images["filter"] == "F277W")].iloc[0]
    f = abs(i.new_minus_old_ap_njy)
    e = i.local_scatter_njy
    z = h.EAZY_z_a
    lo = h.EAZY_l68
    hi = h.EAZY_u68
    m = 31.4 - 2.5 * np.log10(f)
    cs = c[c.site == site]
    fnames = ",".join(sorted(set(cs["filter"])))
    d = flux[(flux.site == site) & (flux["filter"] == "F277W")]
    d1 = d[d.pid == 8060].f15_njy
    d2 = d[d.pid == 12577].f15_njy
    records.append(
        dict(
            site=site,
            **sites.loc[site].to_dict(),
            pointlike_filters=fnames,
            dr2_id=int(h.dr2_id),
            host_sep_arcsec=h.sep_arcsec,
            zphot=z,
            zphot_68lo=lo,
            zphot_68hi=hi,
            difference_f277w_ap_njy=i.new_minus_old_ap_njy,
            diff_local_scatter_njy=e,
            abs_delta_ab_ap=fn(m, z),
            bright_envelope=fn(31.4 - 2.5 * np.log10(f + e), hi),
            faint_envelope=fn(31.4 - 2.5 * np.log10(f - e), lo),
            rest_band_center_um=2.77 / (1 + z),
            old_dither_median_njy=np.median(d1),
            new_dither_median_njy=np.median(d2),
            old_dither_count=len(d1),
            new_dither_count=len(d2)
        )
    )
pd.DataFrame(records).to_csv("nuclear_candidates.csv", index=False)
print(
    pd.DataFrame(records)
    .set_index("site")[
        ["abs_delta_ab_ap", "bright_envelope", "faint_envelope", "rest_band_center_um"]
    ]
    .round(2)
)
