"""Plot closest approaches from the frozen Spitzer and WISE epoch fits."""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_plot(folder, output):
    forecast = load_module("spitzer_forecast", folder / "forecast.py")
    wise = load_module("wise_epoch_fit", folder / "wise-check" / "compare.py")
    lens = pd.read_csv(folder / "lenses.csv").set_index("name").loc[forecast.NAME]
    src = pd.read_csv(
        folder / "source_covariance.csv", dtype={"source_id": str}
    ).set_index("source_id")
    epochs = pd.read_csv(folder / "spitzer_epochs.csv")
    spitzer, _, _ = forecast.fit_spitzer(epochs, lens, np.ones(len(epochs), dtype=bool))
    refs = pd.read_csv(folder / "wise-check" / "reference-offsets.csv")
    w2_epochs = pd.read_csv(folder / "wise-check" / "W2-epochs.csv")
    w2, _ = wise.fit(w2_epochs, refs[refs["product"] == "W2"], True)
    # Convert the WISE fit to the offset/parallax/PM order and J2016 epoch.
    w2_j2016 = np.r_[w2[:2] + (2016 - wise.EPOCH) * w2[2:4], w2[4], w2[2:4]]
    axes = forecast.unit_vectors(lens.ra, lens.dec)
    fig, plots = plt.subplots(1, 2, figsize=(10.2, 4.6), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.22, top=0.78, wspace=0.10)
    models = [
        ("Spitzer", spitzer, "#175781", "-"),
        ("WISE W2", w2_j2016, "#b25027", "--"),
    ]
    for ax, (sid, start, stop), title in zip(
        plots,
        [(forecast.IDS[0], 2033.0, 2035.25), (forecast.IDS[1], 2039.0, 2041.25)],
        ["First background star", "Second background star"],
    ):
        years = forecast.grid(start, stop, 0.5)
        dates = np.datetime64("2000-01-01T12:00") + np.round(
            (years - 2000) * 365.25 * 1440
        ).astype("timedelta64[m]")
        factors = -axes[:2] @ forecast.earth_position(years).T
        state = forecast.params(src.loc[sid], lens, 2016)
        state[2] = max(0.0, state[2])
        for label, fg, color, style in models:
            sep = np.hypot(*forecast.tracks(state - fg, years, factors))
            ax.plot(dates, sep, color=color, ls=style, lw=2, label=label)
            idx = np.argmin(sep)
            ax.scatter(
                dates[idx],
                sep[idx],
                s=30,
                color=color,
                edgecolors="white",
                linewidths=0.7,
                zorder=3,
            )
        ax.set_title(
            title, loc="left", family="DejaVu Sans", fontsize=13, weight="bold", pad=22
        )
        ax.text(
            0,
            1.028,
            "Gaia DR3 " + sid,
            transform=ax.transAxes,
            fontsize=9,
            color="#425163",
        )
        ax.set_xlim(dates[0], dates[-1])
        ax.set_ylim(0, 680)
        ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
        ax.set_yticks([0, 100, 200, 300, 400, 500, 600])
        ax.grid(axis="y", color="#dbe0e4", lw=0.7)
        ax.tick_params(colors="#33414b", labelsize=9)
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#abb5bd")
    plots[0].set_ylabel("Predicted angular separation (mas)", fontsize=11, labelpad=9)
    fig.legend(
        *plots[0].get_legend_handles_labels(),
        loc="lower center",
        bbox_to_anchor=(0.53, 0.015),
        ncol=2,
        frameon=False,
        fontsize=11,
    )
    fig.suptitle(
        "WISE J2000+3629 approaches two background stars",
        x=0.09,
        ha="left",
        y=0.968,
        fontsize=17,
        weight="bold",
    )
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("encounters.png")
    )
    args = parser.parse_args()
    make_plot(args.data, args.output)
