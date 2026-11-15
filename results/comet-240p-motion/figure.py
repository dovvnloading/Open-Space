"""Plot the orbital time lag and departure from the local tangent."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from observations import match_components, read_mpc
from replay import HERE, REF, Motion, solve


def plot():
    motion = Motion(HERE / "data" / "states.csv")
    pairs = match_components(read_mpc(HERE / "data" / "astrometry.mpc"))
    mpc = pairs.groupby(["station", "night"])[["jd", "E", "N"]].median().reset_index()
    nordic = pd.read_csv(HERE / "data" / "not.csv")
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.3), layout="constrained")
    for label, sample, color in [
        ("MPC, 40 sites", mpc, "#20658e"),
        ("Nordic Optical Telescope", nordic, "#b55e23"),
    ]:
        design, _ = motion.design(sample)
        parameters = solve(design[:, :6], sample[["E", "N"]].to_numpy().ravel())
        times = np.linspace(max(sample.jd.min(), 2460959.5), min(sample.jd.max(), 2461128.5), 130) - REF
        evolution = []
        for time in times:
            r, v = (motion.matrix(time)[:6, :6] @ parameters).reshape(2, 3)
            evolution.append(motion.observable(r, v, epoch=time))
        quantities = pd.DataFrame(evolution)
        dates = pd.to_datetime(times + REF, unit="D", origin="julian")
        axes[0].plot(dates, quantities.trailing_minutes, color=color, label=label, lw=2)
        axes[1].plot(dates, quantities.offset_from_local_orbital_tangent_km / 1000, color=color, lw=2)

    ticks = pd.to_datetime(["2025-11-01", "2026-01-01", "2026-03-01"])
    for axis, title, units in zip(axes, ["Time behind component A", "Distance off A’s local orbital tangent"],
                                  ["minutes", "thousand km"]):
        axis.set_title(title, fontsize=11, loc="left")
        axis.set_ylabel(units)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=.25)
        axis.set_xticks(ticks, ["Nov 2025", "Jan 2026", "Mar 2026"])
    axes[0].legend(fontsize=8, frameon=False)
    figure.savefig(HERE / "lag.png", dpi=150)


if __name__ == "__main__":
    plot()
