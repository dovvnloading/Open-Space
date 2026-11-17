"""Visualize which systemic velocities would keep the strongest DESI candidate unbound."""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from infer import draws, escape, load_inputs
from screen_desi import cartesian
from make_queue import number

if __name__ == "__main__":
    stars, gs = load_inputs()
    s = next(s for s in stars if s["source_id"] == "1717410430709663616")
    g = gs[s["source_id"]]
    sky = (number(g["ra"]), number(g["dec"]), 1.0, 0.0, 0.0)
    solar = cartesian(*sky, 0.0)[1]
    unit = cartesian(*sky, 1.0)[1] - solar
    grid = np.linspace(220, 535, 320)
    fig, ax = plt.subplots(figsize=(8.2, 5.0), layout="constrained")
    colors = ["#214e80", "#b85b31"]
    for model, style, label in [
        ("spec_plx", "-", "SpecDis + parallax"),
        ("bj_photogeo", "--", "Bailer–Jones photogeometric"),
    ]:
        _, x, v = draws(g, s, model)
        proj = unit @ v
        vt2 = np.sum(v * v, axis=0) - proj**2
        for e, col, name in zip(escape(x), colors, ["McMillan 2017", "Gala 2022"]):
            p = np.mean(
                vt2[:, None] + (grid + unit @ solar) ** 2 > e[:, None] ** 2, axis=0
            )
            ax.plot(grid, p, style, color=col, lw=2.1, label=f"{name} · {label}")
    ax.axvspan(
        number(s["rv"]) - 4.65, number(s["rv"]) + 4.65, color="#5d8b64", alpha=0.22
    )
    ax.axvline(number(s["rv"]), color="#41724a", lw=1.5)
    ax.text(490, 0.07, "Measured\nepoch RV", ha="right", va="bottom", color="#355c3e")
    ax.set(
        xlim=(220, 535),
        ylim=(-0.02, 1.05),
        xlabel="Assumed systemic heliocentric radial velocity (km/s)",
        ylabel="Fraction above local escape speed",
        title="Gaia DR3 1717410430709663616",
    )
    fig.legend(loc="outside lower center", ncols=2, frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.18)
    fig.savefig("systemic_velocity.png", dpi=160)
