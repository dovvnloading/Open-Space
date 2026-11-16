"""Map cellwise OLCI changes and compare paired daily spectra on their common footprint."""

import csv
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from aaot_sensitivity import DEPTH_URL
from matplotlib import colors
from matplotlib.ticker import FuncFormatter
from olci_zarr import read
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent


def load_table(file):
    with (ROOT / file).open() as f:
        return list(csv.DictReader(f))


def main():
    rows = load_table("aaot_paired_spectra.csv")
    paired = load_table("aaot_spectral_context.csv")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": "#304754",
            "xtick.color": "#304754",
            "ytick.color": "#304754",
            "text.color": "#1e3340",
            "axes.titleweight": "bold",
            "axes.edgecolor": "#8296a1",
            "figure.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(11.5, 5.7))
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1, 1.38],
        left=0.07,
        right=0.98,
        bottom=0.24,
        top=0.76,
        wspace=0.28,
    )
    ax = fig.add_subplot(gs[0, 0])
    bx = fig.add_subplot(gs[0, 1])
    fig.text(
        0.07,
        0.93,
        "A one-day red reflectance peak off Venice",
        fontsize=21,
        fontweight="bold",
    )
    fig.text(
        0.07,
        0.865,
        "Sentinel-3 OLCI  ·  Gulf of Venice  ·  5–7 November 2025",
        fontsize=12,
        color="#4b606c",
    )
    with netcdf_file(io.BytesIO(read(DEPTH_URL)), mmap=False) as etopo:
        lon = etopo.variables["longitude"].data.copy()
        lat = etopo.variables["latitude"].data.copy()
        alt = etopo.variables["altitude"].data.copy()
    ax.set_facecolor("#eef3f5")
    ax.contourf(lon, lat, alt, levels=[0, 4000], colors=["#d2d0ca"])
    ax.contour(lon, lat, alt, levels=[0], colors=["#8d928f"], linewidths=0.55)
    selected = [r for r in rows if float(r["depth_m"]) > 0]
    xs = np.array([float(r["lon"]) for r in selected])
    ys = np.array([float(r["lat"]) for r in selected])
    dz = np.array([float(r["delta_C"]) for r in selected])
    np.unique(xs)
    np.unique(ys)
    allx = np.arange(round((xs.max() - xs.min()) * 24) + 1) / 24 + xs.min()
    ally = np.arange(round((ys.max() - ys.min()) * 24) + 1) / 24 + ys.min()
    z = np.full((len(ally), len(allx)), np.nan)
    z[
        np.rint((ys - ys.min()) * 24).astype(int),
        np.rint((xs - xs.min()) * 24).astype(int),
    ] = dz
    palette = colors.LinearSegmentedColormap.from_list(
        "change", ["#16547f", "#68a3c5", "#f6f4ee", "#d69366", "#a84723"]
    )
    norm = colors.TwoSlopeNorm(vmin=-0.25, vcenter=0, vmax=0.25)
    kw = {"cmap": palette, "norm": norm, "shading": "flat", "rasterized": True}
    mesh = ax.pcolormesh(
        np.append(allx - 1 / 48, allx[-1] + 1 / 48),
        np.append(ally - 1 / 48, ally[-1] + 1 / 48),
        z,
        **kw,
    )
    center = (12.5083, 45.3139)
    angle = np.linspace(0, 2 * np.pi, 360)
    ax.plot(
        center[0] + 50 * np.cos(angle) / (111.195 * np.cos(np.radians(center[1]))),
        center[1] + 50 * np.sin(angle) / 111.195,
        "--",
        linewidth=0.9,
        color="#304754",
    )
    ax.plot(
        *center,
        marker="*",
        ms=12,
        color="#ffdc75",
        markeredgecolor="#1f3440",
        markeredgewidth=0.8,
        zorder=5,
    )
    ax.annotate(
        "Acqua Alta tower",
        center,
        xytext=(12.39, 45.77),
        fontsize=8,
        fontweight="bold",
        ha="center",
        arrowprops={"arrowstyle": "-", "color": "#233b49", "lw": 0.8},
    )
    ax.text(12.84, 44.91, "50 km", color="#304754", fontsize=8, ha="left")
    ax.set_xlim(11.88, 13.67)
    ax.set_ylim(44.64, 45.88)
    ax.set_aspect(1 / np.cos(np.radians(center[1])))
    ax.set_xticks([12, 12.5, 13, 13.5])
    ax.set_yticks([44.75, 45, 45.25, 45.5, 45.75])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}°E"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}°N"))
    ax.tick_params(length=2)
    ax.set_title("Change in 620 nm curvature, 6–7 Nov", loc="left", fontsize=11, pad=9)
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, pos.y0 - 0.092, pos.width, 0.017])
    cb = fig.colorbar(
        mesh,
        cax=cax,
        orientation="horizontal",
        ticks=[-0.25, -0.125, 0, 0.125, 0.25],
        extend="both",
    )
    cb.outline.set_linewidth(0.4)
    cb.ax.set_xlabel("Change in log curvature")
    cb.ax.tick_params(length=2, labelsize=8)
    dates = [
        ("2025-11-05", "2025-11-06", "#bf5d2b", "o", "5–6 Nov"),
        ("2025-11-06", "2025-11-07", "#146185", "s", "6–7 Nov"),
        ("2025-11-05", "2025-11-07", "#536067", "D", "5–7 Nov"),
    ]
    waves = [443, 490, 510, 560, 620, 665, 674, 681, 709]
    for a, b, color, m, label in dates:
        r = next(r for r in paired if r["start"] == a and r["end"] == b)
        vals = [100 * np.expm1(float(r[f"median_log_change_{w}"])) for w in waves]
        bx.plot(
            waves,
            vals,
            label=f"{label}  ({r['paired_all_bands']} cells)",
            color=color,
            marker=m,
            ms=4,
            linestyle="--" if label == "5–7 Nov" else "-",
            lw=1.6,
        )
    bx.axhline(0, color="#88969a", lw=0.8, zorder=0)
    bx.grid(axis="y", color="#dce2e3", lw=0.6)
    bx.set_axisbelow(True)
    bx.set_xticks([443, 490, 560, 620, 665, 709])
    bx.set_xlim(433, 717)
    bx.set_xlabel("OLCI band center (nm)")
    bx.set_ylabel("Median paired reflectance change (%)")
    bx.set_title(
        "Daily reflectance change within 50 km", loc="left", fontsize=11, pad=9
    )
    bx.legend(loc="upper left", frameon=False, fontsize=8)
    fig.text(
        0.07,
        0.03,
        "Map shows wet cells observed on both days. Spectra use cells positive at all nine bands on both dates.\nModified Copernicus Marine Service information (DOI 10.48670/moi-00280). Coastline from NOAA ETOPO1.",
        fontsize=8,
        color="#516371",
    )
    fig.savefig(ROOT / "aaot_optical_pulse.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    main()
