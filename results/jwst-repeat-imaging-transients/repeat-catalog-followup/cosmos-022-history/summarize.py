"""Collapse only distinct UTC nights and filters for calibrated COSMOS-022 data."""

from pathlib import Path
from hashlib import sha256
import numpy as np, pandas as pd, json, sys
from matplotlib import pyplot as plt
import matplotlib.dates as dates
from datetime import datetime


def mad(x):
    return 1.4826 * np.nanmedian(abs(x - np.nanmedian(x))) if len(x) > 1 else np.nan


def main():
    singles = pd.read_csv("phot_l2.csv")
    mosaics = pd.read_csv("phot_l3.csv")
    for level, r in [(2, singles), (3, mosaics)]:
        query = pd.read_csv(f"level{level}.csv")
        assert set(query.obs_id) == set(r.obs_id) and r.obs_id.is_unique
    g = singles.query('status=="measured"').copy()
    g["night"] = g.utc.str[:10]
    rows = []
    for (night, band), d in g.groupby(["night", "filter"]):
        rows.append(
            dict(
                utc_night=night,
                filter=band,
                pids=",".join(sorted(set(d.obs_id.str[:7]))),
                n=len(d),
                mjd_first=d.mjd_start.min(),
                mjd_last=d.mjd_end.max(),
                median_fixed_njy=d.fixed_njy.median(),
                dither_mad_njy=mad(d.fixed_njy),
                min_njy=d.fixed_njy.min(),
                max_njy=d.fixed_njy.max(),
                median_local_blank_scatter_njy=d.local_scatter_njy.median(),
                median_centered_njy=d.centered_njy.median(),
                median_wide_njy=d.wide_njy.median(),
                median_east_arcsec=d.east_arcsec.median(),
                median_north_arcsec=d.north_arcsec.median(),
            )
        )
    tab = pd.DataFrame(rows)
    tab.to_csv("nightly.csv", index=False, float_format="%.9g")
    m = mosaics.query('status=="measured"').copy()
    m["month"] = m.utc.str[:7]
    pivot = m.groupby(["filter", "month"]).fixed_njy.median().unstack()
    v = pivot.loc[:, ["2025-12", "2026-05"]].dropna()
    v.columns = ["old_njy", "new_njy"]
    v["delta_njy"] = v.new_njy - v.old_njy
    v["delta_ab_mag"] = -2.5 * np.log10(v.new_njy / v.old_njy)
    v.to_csv("medium_changes.csv", float_format="%.8g")
    p = []
    for n in ["level2.csv", "level3.csv", "mast_cone.csv"]:
        b = Path(n).read_bytes()
        p.append(dict(input=n, sha256=sha256(b).hexdigest(), bytes=len(b)))
    pd.DataFrame(p).to_csv("inventory_checksums.csv", index=False)
    groups = ["F150W", "F200W", "F140M", "F162M", "F182M", "F210M"]
    colors = ("#625b71", "#302e36", "#d25922", "#973755", "#1b718b", "#3d7b42")
    markers = ["o", "s", "P", "D", "^", "v"]
    fig, ax = plt.subplots(
        1, 2, figsize=(9.2, 3.8), width_ratios=[1.35, 1], gridspec_kw={"wspace": 0.24}
    )
    for a, mindate in zip(ax, ["2022-01-01", "2025-12-01"]):
        for fil, col, marker in zip(groups, colors, markers):
            t = tab[(tab["filter"] == fil) & (tab.utc_night >= mindate)]
            t0 = pd.to_datetime(t.utc_night)
            a.errorbar(
                t0,
                t.median_fixed_njy,
                yerr=t.dither_mad_njy.fillna(0),
                fmt=marker,
                color=col,
                ms=4,
                lw=0.8,
                capsize=1.5,
                label=fil,
            )
        a.grid(axis="y", alpha=0.22)
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=8)
        a.set_ylim(
            tab[tab["filter"].isin(groups)].median_fixed_njy.min() - 55,
            tab[tab["filter"].isin(groups)].median_fixed_njy.max() + 60,
        )
    ax[0].xaxis.set_major_locator(dates.YearLocator())
    ax[0].xaxis.set_major_formatter(dates.DateFormatter("%Y"))
    ax[0].set_ylabel("Host-inclusive 0.16″ aperture (nJy)")
    ax[0].set_title("Multi-year near-infrared flux", loc="left", fontsize=11)
    ax[1].xaxis.set_major_locator(dates.MonthLocator(interval=2))
    ax[1].xaxis.set_major_formatter(dates.DateFormatter("%b %y"))
    ax[1].set_title("2025–26 medium-band visits", loc="left", fontsize=11)
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.20, top=0.83)
    fig.suptitle("COSMOS-022", x=0.075, ha="left", fontweight="bold", fontsize=14)
    fig.savefig("cosmos022-lightcurve.png", dpi=155)
    plt.close(fig)


if __name__ == "__main__":
    main()
