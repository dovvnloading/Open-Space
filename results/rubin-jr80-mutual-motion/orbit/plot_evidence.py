"""Plot the Rubin and HST relative astrometry and an independent-epoch fit check."""

import csv
from datetime import date

import matplotlib.pyplot as plt
import numpy as np

from fit import BASE, PARAMETERS, predict, rows


def records(path):
    with path.open() as file:
        return list(csv.DictReader(file))


def axis_line(ax, east, north, color):
    ax.plot([0, east], [0, north], color=color, lw=1.3, alpha=0.65)
    ax.scatter(east, north, color=color, s=42, zorder=4, edgecolors="white", lw=0.6)


def main():
    root = BASE.parent
    fig, (sky, check) = plt.subplots(
        1, 2, figsize=(10.2, 4.9), gridspec_kw={"width_ratios": [1, 1.12]}
    )
    fig.subplots_adjust(wspace=0.36, left=0.09, right=0.98, top=0.85, bottom=0.18)
    palette = ["#315a82", "#2b876f", "#b56533", "#934063"]
    hst = records(root / "hst" / "measurements.csv")
    e = np.mean([float(r["core_east_arcsec"]) for r in hst])
    n = np.mean([float(r["core_north_arcsec"]) for r in hst])
    axis_line(sky, e, n, palette[0])
    sky.text(e + 0.02, n, "June 13 · HST", va="center", fontsize=9)
    raw = records(root / "output" / "astrometry.csv")
    for color, r in zip(palette[1:], (r for r in raw if r["band"] == "i")):
        e, n = float(r["east"]), float(r["north"])
        axis_line(sky, e, n, color)
        date_text = (
            date.fromisoformat(r["utc_date"]).strftime("%b %d").replace(" 0", " ")
        )
        sky.text(e + 0.025, n, date_text, va="center", fontsize=9)
    sky.scatter(0, 0, marker="+", s=70, color="#303941")
    sky.set_aspect("equal")
    sky.set(
        xlim=(-0.08, 0.88),
        ylim=(-0.25, 0.95),
        xlabel="East offset (arcsec)",
        ylabel="North offset (arcsec)",
    )
    sky.set_title("Relative positions in 2026", loc="left", fontsize=12, weight="bold")
    # These fits use only archival CFHT and HST; the three Rubin points are withheld.
    fits = records(BASE / "fit_no_rubin.csv")
    measured = np.array(
        [[float(r["east_arcsec"]), float(r["north_arcsec"])] for r in rows]
    )
    include = np.array([r["id"].startswith("2026-") for r in rows])
    if sum(include) != 3:
        raise ValueError("Expected three Rubin observing nights")
    counter = np.arange(3)
    for f, color, marker in zip(
        fits[:3], ["#1e7766", "#b56434", "#797d90"], ["o", "s", "^"]
    ):
        p = [float(f[k]) for k in PARAMETERS]
        predicted = predict(p)[include]
        obs = measured[include]
        error = np.minimum(
            np.linalg.norm(predicted - obs, axis=1),
            np.linalg.norm(predicted + obs, axis=1),
        )
        check.plot(
            counter,
            error,
            marker=marker,
            color=color,
            lw=1.7,
            ms=6,
            label=f"{float(f['period_days']):.1f} days",
        )
    check.set(
        xticks=counter,
        xticklabels=["Jun 27", "Jul 7", "Jul 9"],
        ylim=(0, 0.17),
        ylabel="Relative-position error (arcsec)",
        xlabel="Rubin observation (UTC)",
    )
    check.set_title(
        "Predicting the withheld Rubin nights", loc="left", fontsize=12, weight="bold"
    )
    check.legend(
        title="CFHT + HST period",
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    for ax in (sky, check):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#d8dee2", alpha=0.52, lw=0.6)
    fig.suptitle("(556068) 2014 JR80", x=0.09, ha="left", size=18, weight="bold")
    fig.savefig(root / "output/orbit_evidence.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
