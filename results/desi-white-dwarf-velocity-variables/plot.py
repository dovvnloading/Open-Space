"""Plot selected exposure velocities. Points at each exposure have small display offsets."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
exposures = pd.read_csv(ROOT / "selected-exposures.csv")
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.5), sharey=True)
for ax, name in zip(axes, ["WDJ023809.44+273329.79", "WDJ042155.26-132531.35"]):
    table = exposures[exposures.WDJname == name].sort_values("BMJD")
    offset = int(name.startswith("WDJ0238"))
    x = np.arange(len(table)) + offset
    labels = [str(int(value)) for value in table.EXPOSURE]
    if offset:
        labels.insert(0, "LAMOST\n2014")
        header = fits.getheader(ROOT / "independent" / "lamost-257903093.fits.gz")
        ax.errorbar(0, 299792.458*header["Z"], yerr=299792.458*header["Z_ERR"],
                    fmt="D", color="#4b4b4b", capsize=3)
    for arm, displacement, color, symbol in [("B", -.06, "#155c99", "o"), ("R", .06, "#bc482d", "s")]:
        ax.errorbar(x + displacement, table["RV_"+arm], yerr=table["RV_"+arm+"_error"],
                    fmt=symbol, color=color, capsize=3, label={"B":"DESI blue", "R":"DESI red"}[arm])
    ax.axhline(0, color="#8c9499", linewidth=.8)
    ax.set_title(name.replace("WDJ","WD J"), fontsize=12)
    ax.set_xticks(np.arange(len(labels)), labels, fontsize=9)
    ax.set_xlabel("Exposure in chronological order")
    ax.grid(axis="y", alpha=.2)
    ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("Apparent radial velocity (km/s)")
axes[1].legend(frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(ROOT / "velocity-states.png", dpi=160, bbox_inches="tight")
