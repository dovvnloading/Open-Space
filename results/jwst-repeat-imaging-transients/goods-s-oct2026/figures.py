"""Aperture light curves and matched-sky-pixel comparisons for the reviewed sources."""

from pathlib import Path
from datetime import datetime
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.dates as dates
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from scipy.ndimage import map_coordinates
import warnings, numpy as np, pandas as pd

warnings.simplefilter("ignore", FITSFixedWarning)
plt.rcParams.update(
    {
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)
pos = pd.concat(
    [pd.read_csv("sites.csv"), pd.read_csv("onsets.csv")[["site", "ra", "dec"]]]
).set_index("site")
sh = pd.read_csv("catalog_stats.csv").set_index(["obs", "filter"])
dy, dx = np.mgrid[-25:26, -25:26] * 0.03


def get(site, pid, band):
    obs = (
        int(pd.read_csv("multiband.csv").query("group==@site").obs.iloc[0])
        if site < 1000
        else 1
    )
    root = "stamps" if site < 1000 else "follow_stamps"
    z = np.load(f"{root}/{pid}-{obs}-{band}.npz")
    j = list(z["site"]).index(site)
    h = fits.Header.fromstring(str(z["header"]), sep="\n")
    w = WCS(h, naxis=2)
    ra, dec = pos.loc[site, ["ra", "dec"]]
    d = sh.loc[obs, f"F{band}W"]
    east, north = (
        (d.offset_east_arcsec, d.offset_north_arcsec) if pid == 8060 else (0, 0)
    )
    x, y = w.all_world2pix(
        ra + (dx - east) / 3600 / np.cos(np.deg2rad(dec)), dec + (dy - north) / 3600, 0
    )
    return map_coordinates(
        z["data"][j],
        [y - z["y0"][j], x - z["x0"][j]],
        order=1,
        mode="constant",
        cval=np.nan,
    )


def display(ax, im, lo, hi, title, side=0.765):
    ax.imshow(
        im,
        origin="lower",
        extent=[side, -side, -side, side],
        cmap="gray",
        vmin=lo,
        vmax=hi,
    )
    ax.set_title(title, fontsize=9.5)
    ax.set(xticks=[], yticks=[])
    ax.scatter(0, -0.19, marker="^", color="#fa845a", s=13)
    ax.text(
        0.02,
        0.04,
        "N ↑  E ←",
        transform=ax.transAxes,
        fontsize=6.5,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.45, "pad": 1, "edgecolor": "none"},
    )


def light(ax, site, df):
    bands = ["115", "150", "200", "277", "356", "444"]
    colors = ["#54596e", "#307bb9", "#21947b", "#d06531", "#743482", "#9f2549"]
    marks = ["d", "^", "v", "o", "s", "P"]
    q = df[df.site == site]
    for b, col, m in zip(bands, colors, marks):
        a = q[q["filter"] == f"F{b}W"].sort_values("utc")
        ax.errorbar(
            pd.to_datetime(a.utc),
            a.aperture_njy,
            yerr=a.local_scatter_njy,
            fmt=m,
            ms=4,
            elinewidth=1,
            capsize=1.3,
            color=col,
            label=f"F{b}W",
        )
    ax.xaxis.set_major_locator(dates.YearLocator())
    ax.xaxis.set_major_formatter(dates.DateFormatter("%Y"))
    ax.grid(axis="y", alpha=0.2)
    ax.tick_params(axis="x", rotation=23, labelsize=7.5)
    ax.set_ylabel("")
    ax.axhline(0, color="gray", lw=0.6)


for site in [185, 99, 1004]:
    fig = plt.figure(figsize=(9.5, 5.2 if site == 1004 else 4.65))
    gs = fig.add_gridspec(
        2,
        4,
        width_ratios=[1, 1, 1, 1.8],
        height_ratios=[1, 1],
        wspace=0.24,
        hspace=0.28,
    )
    axlc = fig.add_subplot(gs[:, 3])
    df = pd.read_csv(
        "follow_forced_lightcurves.csv" if site == 1004 else "forced_lightcurves.csv"
    )
    df = pd.concat([df, pd.read_csv("history_forced_lightcurves.csv")])
    light(axlc, site, df)
    axlc.set_title("0.15″ aperture flux (nJy)")
    axlc.legend(
        ncol=2,
        fontsize=7,
        loc="best",
        handlelength=0.9,
        columnspacing=0.45,
        frameon=False,
    )
    for i, b in enumerate(["277", "356"] if site != 1004 else ["277"]):
        a = get(site, 8060, b)
        n = get(site, 12577, b)
        d = n - a
        center = np.hypot(dx, dy) < 0.3
        lo = np.nanquantile(np.r_[a.flatten(), n.flatten()], 0.20)
        hi = max(np.nanmax(n[center]), np.nanmax(a[center]))
        dd = np.nanquantile(np.abs(d[center]), 0.99)
        plots = [
            (a, lo, hi, "October 2025"),
            (n, lo, hi, "October 2026"),
            (d, -dd, dd, "2026 minus 2025"),
        ]
        for k, (im, v1, v2, label) in enumerate(plots):
            display(fig.add_subplot(gs[i, k]), im[:, ::-1], v1, v2, f"F{b}W · {label}")
    if site == 1004:
        ff = sorted(
            (
                p
                for p in Path("single_stamps").glob("*.npz")
                if str(np.load(p)["filter"]) == "F277W"
            ),
            key=lambda p: str(np.load(p)["date"]),
        )
        for k, p in enumerate(ff):
            z = np.load(p)
            display(
                fig.add_subplot(gs[1, k]),
                z["im"][::-1, ::-1][::-1, :],
                -0.01,
                0.35,
                f"Dither {k+1} · " + str(z["date"])[11:16],
                1.26,
            )
    fig.suptitle(
        f"GOODS-S {site}"
        + ("  |  Fixed detector pixel" if site == 1004 else "  |  Nuclear variability"),
        x=0.04,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.04,
        0.025,
        "Orange marker indicates the tested sky position. Stamps in each epoch pair use the same surface-brightness scale.",
        fontsize=8,
    )
    fig.subplots_adjust(top=0.88, bottom=0.15, left=0.035, right=0.985)
    fig.savefig(f"goods-{site}.png", dpi=165)
    plt.close(fig)
