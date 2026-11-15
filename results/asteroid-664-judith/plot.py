"""Draw the diagnostic from the committed per-transit predictions."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
data = pd.read_csv(root / "measurements/transits.csv")
selected = data.loc[data.catalog_window].copy()
elapsed = (selected.epoch - selected.epoch.min()).to_numpy() * 24
x = (elapsed - elapsed.mean()) / 24
basis = np.column_stack([np.ones(len(x)), x, x**2])
def remove_trend(s):
    return s - basis @ np.linalg.lstsq(basis, s, rcond=None)[0]

fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9), sharey=True, constrained_layout=True)
for ax, model, color in zip(axes, ("5664", "5665"), ("#10684e", "#a64518")):
    prediction = remove_trend(selected[f"al_{model}"].to_numpy())
    observed = remove_trend(selected.fpr_oc.to_numpy())
    ax.scatter(prediction, observed, s=34, c=color)
    ax.axhline(0, color="#d0d0d0", lw=0.8)
    ax.axvline(0, color="#d0d0d0", lw=0.8)
    ax.plot([-1.55, 1.55], [-1.55, 1.55], color="#686868", ls="--", label="Equal displacement")
    ax.set(xlim=(-1.55, 1.55), ylim=(-1.55, 1.55), xlabel="Predicted displacement (mas)",
           title=f"DAMIT {model}   ·   r = {np.corrcoef(prediction, observed)[0, 1]:+.3f}")
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Gaia displacement (mas)")
axes[0].legend(frameon=False, loc="upper left", fontsize=8)
fig.suptitle("(664) Judith  |  Gaia positions, 13–16 December 2017", fontsize=13)
fig.savefig(root / "measurements/pole-comparison.png", dpi=150)
