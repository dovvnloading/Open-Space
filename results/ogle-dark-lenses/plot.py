"""Plot an OGLE event with its no-parallax and parallax fits.

Usage: python plot.py DATADIR EVENT OUT.png [--tmax HJD-2450000] [--forecast DAYS]
"""
import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import candidate
import lens


def binned(t, y, width):
    edges = np.arange(t.min(), t.max() + width, width)
    index = np.digitize(t, edges)
    out = [(t[index == i].mean(), np.median(y[index == i]), 1.4826 * np.median(np.abs(y[index == i] - np.median(y[index == i]))) / np.sqrt((index == i).sum()))
           for i in np.unique(index) if (index == i).sum() >= 3]
    return np.array(out).T


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir")
    parser.add_argument("event")
    parser.add_argument("out")
    parser.add_argument("--tmax", type=float, default=None)
    parser.add_argument("--forecast", type=float, default=0.0)
    args = parser.parse_args()
    data = candidate.prepare(pathlib.Path(args.datadir), args.event, args.tmax)
    t, flux, sigma, par = data["t"], data["flux"], data["sigma"], data["par"]
    grid = np.arange(t.min(), t.max() + args.forecast, 1.0)
    grid_par = lens.Parallax(grid, data["ra"], data["dec"], data["t0par"])
    mag = lens.ZP - 2.5 * np.log10(np.maximum(flux, 1e-6))
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=False, gridspec_kw={"height_ratios": [3, 1.2, 1.2]})
    axes[0].errorbar(t, mag, yerr=1.0857 * sigma / np.maximum(flux, 1e-6), fmt=".", ms=2, color="0.6", elinewidth=0.3,
                     zorder=1)
    curves = {}
    base = data["base"][0]
    amp = lens.magnification(sum(x * x for x in lens.trajectory(base, grid)))
    _, fs, fb = lens.model(base, t, flux, sigma)
    curves["no parallax"] = (fs * amp + fb, lens.model(base, t, flux, sigma)[0], "C0")
    for sign, color in ((1, "C3"), (-1, "C2")):
        theta = data["fits"][sign][0]
        _, fs, fb = lens.model(theta, t, flux, sigma, par)
        amp = lens.magnification(sum(x * x for x in lens.trajectory(theta, grid, grid_par)))
        label = f"parallax u0{'>' if sign > 0 else '<'}0  tE={theta[2]:.0f} d  piE=({theta[3]:+.3f},{theta[4]:+.3f})"
        curves[label] = (fs * amp + fb, lens.model(theta, t, flux, sigma, par)[0], color)
    for label, (curve, _, color) in curves.items():
        axes[0].plot(grid, lens.ZP - 2.5 * np.log10(np.maximum(curve, 1e-6)), color=color, lw=1.2, label=label)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("I (mag)")
    axes[0].legend(fontsize=8)
    axes[0].set_title(f"OGLE-{args.event}   chi2: no parallax {data['base'][2]:.1f}, "
                      f"u0>0 {data['fits'][1][2]:.1f}, u0<0 {data['fits'][-1][2]:.1f}   N={len(t)}")
    zoom = t > data["t0par"] - 700
    for axis, (label, (_, at_data, color)) in zip(axes[1:], list(curves.items())[:2]):
        resid = (flux - at_data) / np.maximum(at_data, 1e-6) * 1.0857
        axis.plot(t[zoom], resid[zoom], ".", ms=2, color="0.75")
        if zoom.sum() > 20:
            bt, by, be = binned(t[zoom], resid[zoom], 10.0)
            axis.errorbar(bt, by, yerr=be, fmt="o", ms=3, color=color)
        axis.axhline(0, color="k", lw=0.6)
        axis.set_ylim(-0.08, 0.08)
        axis.set_ylabel(f"data - model (mag)\n{label.split('  ')[0]}", fontsize=8)
    axes[2].set_xlabel("HJD - 2450000")
    fig.tight_layout()
    fig.savefig(args.out, dpi=110)


if __name__ == "__main__":
    main()
