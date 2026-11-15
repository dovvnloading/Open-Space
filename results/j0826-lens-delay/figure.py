"""Make the two-panel research figure from the frozen reduced tables."""
from replay import ROOT, TARGET, read, curves, np
import matplotlib.pyplot as plt

g = curves(read("gaia-nights.csv"), "bjd_tcb_minus_2400000p5")
z = curves([r for r in read("ztf-nights.csv") if r["system"] == TARGET], "mjd")
fig, ax = plt.subplots(2, 1, figsize=(8, 5.7), layout="constrained")
colors = {"A": "#aa3b36", "B": "#19638b"}
for im in "BA":
    t, m, e = g[im, "g"]
    shift = 210 if im == "B" else 0
    # A constant magnitude offset is solely for plotting. All scores use unshifted magnitudes.
    offset = -1.94 if im == "B" else 0
    ax[0].errorbar(t + shift, m + offset, yerr=e, fmt="o", ms=4, lw=.8,
                   alpha=.83, color=colors[im],
                   label=("B shifted +210 d, −1.94 mag" if im=="B" else "A as observed"))
    t, m, e = z[im, "zg"]
    forplot = np.unique((t // 40).astype(int))
    groups = [m[(t // 40) == q] for q in forplot]
    centers = [np.median(t[(t // 40) == q]) for q in forplot]
    ok = np.array([len(v) >= 3 for v in groups])
    # Subtract each image's initial zero point to show its change, preserving observing times.
    baseline = np.median(m[(t >= 58900) & (t < 59500)])
    ax[1].scatter(np.array(centers)[ok], np.array([np.median(v) for v in groups])[ok] - baseline,
                  color=colors[im], s=16, label=im, alpha=.85)
for axis, title, xlabel, ylabel in zip(ax, ["Gaia DR3 · the same brightening in two images",
                                             "ZTF g · later change in the image flux ratio"],
                                      ["Image-A epoch (BJD TCB − 2,400,000.5)", "Observed MJD"],
                                      ["Gaia G (mag, offset)", "Change from early median (mag)"]):
    axis.invert_yaxis()
    axis.set_title(title, loc="left", fontsize=12)
    axis.set_xlabel(xlabel); axis.set_ylabel(ylabel); axis.legend(frameon=False, fontsize=8)
    axis.grid(alpha=.15); axis.spines[["right","top"]].set_visible(False)
fig.savefig(ROOT / "echo.svg", metadata={"Date": None})
