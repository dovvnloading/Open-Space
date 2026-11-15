"""Render the two June 2003 composition episodes after running replay.py."""
import json
from datetime import timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import LogLocator, FuncFormatter



def main():
    here = Path(__file__).resolve().parent
    measurements = json.loads((here / "measurements.json").read_text())
    frame = pd.read_csv(here / "june_swics.csv", parse_dates=["utc"]).sort_values("utc")
    hourly = pd.read_csv(here / "june_plasma.csv", parse_dates=["utc"])
    hourly["he_percent"] = 100 * hourly.alpha_density / hourly.density
    frame["C5"] = frame.C_den / frame.C6_C5.where(frame.C6_C5 > 0)
    fig, axes = plt.subplots(3, 2, figsize=(9.1, 6.7), sharey="row", sharex="col")
    fig.subplots_adjust(left=.117, right=.975, bottom=.12, top=.82, hspace=.13, wspace=.12)
    fig.suptitle("Fully stripped carbon falls while helium remains abundant", x=.117,
                 ha="left", y=.969, fontsize=16, weight="bold", color="#232e35")
    fig.text(.117, .912, "Ulysses  ·  June 2003  ·  Shading marks the selected three-hour integrations",
             fontsize=10, color="#53626c")
    blue, orange, teal = "#20568a", "#b75419", "#137664"
    for j, episode in enumerate(measurements["episodes"]):
        first = pd.Timestamp(episode["first_utc"]).floor("3h")
        end = pd.Timestamp(episode["last_utc"]).floor("3h") + pd.Timedelta(hours=3)
        start = first - pd.Timedelta(hours=18)
        stop = end + pd.Timedelta(hours=18)
        part = frame[frame.utc.between(start, stop)].copy()
        mid = part.utc.dt.floor("3h") + pd.Timedelta(minutes=90)
        p = hourly[hourly.utc.between(start, stop)]
        for ax in axes[:, j]:
            ax.set_axisbelow(True)
            ax.grid(axis="y", color="#d6dce0", linewidth=.65)
            ax.axvspan(first, end, color="#e4edee", zorder=0)
            ax.spines[["top", "right"]].set_visible(False)
            ax.spines[["left", "bottom"]].set_color("#9ca9ad")
            ax.tick_params(labelsize=8, colors="#364850")
            ax.set_xlim(start, stop)
        kw = {"markersize": 3.3, "linewidth": 1.15}
        axes[0, j].plot(mid, part.C6_C5.where(part.C6_C5 > 0), "o-", color=blue, label="C⁶⁺/C⁵⁺", **kw)
        axes[0, j].plot(mid, part.O7_O6.where(part.O7_O6 > 0), "s-", color=orange, label="O⁷⁺/O⁶⁺", **kw)
        axes[1, j].plot(mid, part.C_den.where(part.C_den > 0), "o-", color=blue, label="C⁶⁺/O⁶⁺", **kw)
        axes[1, j].plot(mid, part.C5.where(part.C5 > 0), "^-", color=teal, label="C⁵⁺/O⁶⁺, inferred", **kw)
        axes[2, j].plot(p.utc, p.he_percent, color=orange, lw=1.3)
        axes[2, j].axhline(100 * measurements["alpha_proton_controls"]["high_charge_fast_control_median"],
                           color="#57666e", linestyle="--", lw=1, label="Control median")
        for ax in axes[:2, j]:
            ax.set_yscale("log")
            ax.set_ylim(.008, 16)
            ax.yaxis.set_major_locator(LogLocator(base=10))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
        axes[2, j].set_ylim(0, 22)
        axes[2, j].xaxis.set_major_locator(mdates.HourLocator(byhour=[0], tz=timezone.utc))
        axes[2, j].xaxis.set_major_formatter(mdates.DateFormatter("%b %d", tz=timezone.utc))
        axes[0, j].set_title(f"{first:%B %d}   ·   {episode['distance_AU']:.2f} AU",
                            loc="left", fontsize=11, weight="bold", pad=8)
    for ax, label in zip(axes[:, 0], ["Charge-state ratio", "Carbon / O⁶⁺", "Alpha / proton (%)"]):
        ax.set_ylabel(label, fontsize=10)
    for i, ax in enumerate(axes[:, 1]):
        ax.legend(loc="upper right", frameon=True, framealpha=.92, edgecolor="none",
                  fontsize=7.6, ncol=2, handlelength=1.6, columnspacing=.8)
    fig.text(.117, .035, "UTC  ·  Carbon and oxygen from SWICS  ·  Helium and protons from SWOOPS",
             color="#53626c", fontsize=8.5)
    fig.savefig(here / "composition.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
