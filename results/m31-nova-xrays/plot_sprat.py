"""Plot the three pipeline spectra and the resolved H-alpha profile."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
data = np.genfromtxt(HERE / "sprat-spectra.csv", names=True, delimiter=",")
wavelength = data["wavelength_angstrom"]
fig, (full, line) = plt.subplots(1, 2, figsize=(9, 3.8), width_ratios=(2.1, 1), layout="constrained")
colors = ("#24648b", "#a34c22", "#157563")
exposures = ("03:48:53", "03:59:02", "04:09:10")
for i, (color, stamp) in enumerate(zip(colors, exposures), 1):
    counts = data[f"spec_ss_{i}"]
    full.plot(wavelength, counts / np.median(counts[(wavelength > 5400) & (wavelength < 5700)]),
              color=color, lw=1, alpha=0.85, label=stamp)
    window = np.abs(wavelength - 6562.8) < 75
    continuum = np.median(counts[window & (np.abs(wavelength - 6562.8) > 46)])
    profile = counts[window] - continuum
    line.plot(wavelength[window], profile / profile.max(), ".-", color=color, lw=1.3, ms=4)
for wave, label in [(4861.33, "Hβ"), (6562.8, "Hα")]:
    full.axvline(wave, c="#767676", ls=":", lw=0.8)
    full.text(wave+20, 0.96, label, transform=full.get_xaxis_transform(), va="top")
full.set(xlim=(4050, 7900), xlabel="Wavelength (Å)", ylabel="Counts / median continuum",
         title="M31N 2017-06c · SPRAT")
full.legend(frameon=False, fontsize=8, title="June 13, 2017 UTC", title_fontsize=8)
line.axvline(6562.8, c="#767676", ls=":", lw=0.8)
line.set(xlabel="Wavelength (Å)", ylabel="Fraction of line peak", title="Hα per exposure")
for ax in (full, line):
    ax.spines[["top", "right"]].set_visible(False)
fig.savefig(HERE / "sprat-profiles.png", dpi=150)
