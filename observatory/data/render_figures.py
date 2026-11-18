"""Replot the MAVEN and STIX dossier figures from HEAD. Requires Matplotlib."""

import csv
from datetime import datetime, timezone
import hashlib
import io
from pathlib import Path
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

DATA = Path(__file__).resolve().parent
ROOT = DATA.parents[1]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none",
    "svg.hashsalt": "Open-Space", "figure.facecolor": "#13151a", "axes.facecolor": "#13151a",
    "text.color": "#f3f1ed", "axes.labelcolor": "#dadce3", "xtick.color": "#a9adb8",
    "ytick.color": "#a9adb8", "axes.edgecolor": "#414650", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "axes.axisbelow": True,
    "grid.color": "#30333c", "grid.linewidth": 0.65,
})


def rows(path):
    raw = subprocess.check_output(["git", "-C", str(ROOT), "show", f"HEAD:{path}"])
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8")))), hashlib.sha256(raw).hexdigest()


def finish(fig, ax, name, title, digest):
    ax.tick_params(length=0, pad=6)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=14, pad=13)
    fig.savefig(DATA / "figures" / f"{name}.svg",
                metadata={"Date": None, "Description": "Source CSV SHA-256 " + digest})
    plt.close(fig)


def mars():
    data, digest = rows("results/maven-2025-sep-events/daily.csv")
    august = [row for row in data if row["date"].startswith("2025-08-")]
    fig, ax = plt.subplots(figsize=(7.6, 3.8), layout="constrained")
    for column, label, color in [("electron", "Electrons", "#f6a772"), ("ion", "Ions", "#88c8e4")]:
        days = [int(row["date"][-2:]) for row in august]
        flux = [float(row["max_hourly_median_" + column]) for row in august]
        ax.scatter(days, flux, label=label, color=color, s=18, zorder=3)
    ax.set(yscale="log", xlabel="Day of August 2025 (UTC)",
           ylabel="Peak hourly median flux [(cm² s sr)⁻¹]", xticks=[1, 5, 10, 15, 20, 25, 31])
    ax.legend(loc="upper left", frameon=False)
    finish(fig, ax, "mars", "Energetic particles measured in Mars orbit", digest)


def sun():
    data, digest = rows("results/2021-july-13-farside-eruption/inputs/stix_rates.csv")
    selected = [row for row in data if row["product"] == "SCI" and row["band_keV"] == "15-25"]
    times = [datetime.fromisoformat(row["spacecraft_utc"]).replace(tzinfo=timezone.utc) for row in selected]
    rates = [float(row["count_rate_s"]) for row in selected]
    fig, ax = plt.subplots(figsize=(7.6, 3.8), layout="constrained")
    ax.plot(times, rates, color="#f4b77a", marker=".", markersize=3, linewidth=1.1)
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 15, 30, 45]))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M", tz=timezone.utc))
    ax.set(xlabel="13 July 2021 (spacecraft UTC)", ylabel="Count rate [counts s⁻¹]", ylim=(0, None))
    finish(fig, ax, "sun", "Hard X-rays at 15–25 keV", digest)


if __name__ == "__main__":
    (DATA / "figures").mkdir(exist_ok=True)
    mars()
    sun()
