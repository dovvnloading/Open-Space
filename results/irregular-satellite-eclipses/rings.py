"""Predict the fraction of direct sunlight reaching Phoebe through Saturn's rings."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from scipy.stats import qmc

import search

PROFILE = (
    "https://pds-rings.seti.org/holdings/volumes/VG_28xx/VG_2801/EASYDATA/KM050/PS1P01."
)
VIMS_BASE = (
    "https://opus.pds-rings.seti.org/holdings/volumes/COVIMS_8xxx/COVIMS_8001/data/"
)
VIMS_NAMES = ["VIMS_2006_252_ALPTAU_I_TAU_10KM", "VIMS_2008_160_GAMCRU_I_TAU_10KM"]
MID = 4444.64364601  # Days TDB since 2020-01-01, approximate geometric mid-event.


def disk(samples):
    radius, phi = qmc.Sobol(2, scramble=False).random_base2(samples).T
    return np.array(
        [
            np.sqrt(radius) * np.cos(2 * np.pi * phi),
            np.sqrt(radius) * np.sin(2 * np.pi * phi),
        ]
    )


def pole_at(t):
    centuries = (search.J0 + t - 2451545) / 36525
    ra, dec = np.deg2rad([40.589 - 0.036 * centuries, 83.537 - 0.004 * centuries])
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def sunlight(states, t, uv, profiles):
    delta, axis, dp, ds, dot = search.shadow(states, 609, 699, t)
    sun = axis * dot + delta
    pole = pole_at(t)
    north = search.unit(pole - sun * np.dot(pole, sun))
    east = np.cross(north, sun)
    rays = sun[:, None] + (east[:, None] * uv[0] + north[:, None] * uv[1]) * np.tan(
        np.arcsin(search.SUN_R / ds)
    )
    # Ray from the moon. Planet center at P=dp*axis. Ring plane through P.
    # Offset at a plane close to the planet avoids subtracting enormous squared distances.
    p = dp * axis
    start = rays * (dp / (axis @ rays)) - p[:, None]
    a, b, *_ = search.PLANETS[699]
    stretch = a * a / (b * b) - 1.0
    sp = pole @ start
    vp = pole @ rays
    aa = np.sum(rays * rays, axis=0) + stretch * vp * vp
    ab = np.sum(start * rays, axis=0) + stretch * sp * vp
    cc = np.sum(start * start, axis=0) + stretch * sp * sp
    clear = cc - ab * ab / aa > a * a
    plane_dist = np.dot(p, pole) / vp
    intercept = rays * plane_dist - p[:, None]
    radial = np.sqrt(np.sum(intercept * intercept, axis=0))
    mu = np.abs(vp) / np.sqrt(np.sum(rays * rays, axis=0))
    weights = 1 - 0.6 * (1 - np.sqrt(np.maximum(0, 1 - uv[0] ** 2 - uv[1] ** 2)))
    values = [np.mean(clear)]
    for i, (first, spacing, opacities) in enumerate(profiles):
        # Thin, azimuthally symmetric ring. Negative noise is clipped to vacuum.
        index = np.rint((radial - first) / spacing).astype(int)
        inside = (index >= 0) & (index < len(opacities)) & (plane_dist > 0)
        tau = np.where(inside, opacities[np.clip(index, 0, len(opacities) - 1)], 0)
        transmission = clear * np.exp(-tau / mu)
        values.append(np.mean(transmission))
        if i == 0:
            values.append(np.average(transmission, weights=weights))
    return values


def run(cache, outdir, samples):
    files = {}
    profiles = []
    for name, base, spacing in [("PS1P01", PROFILE[: -len("PS1P01.")], 50)] + [
        (n, VIMS_BASE, 10) for n in VIMS_NAMES
    ]:
        for ext in ["TAB", "LBL"]:
            fname = f"{name}.{ext}"
            url = base + fname
            p = cache / (fname if name == "PS1P01" else fname.lower())
            if not p.exists():
                p.write_bytes(
                    urlopen(
                        Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=90
                    ).read()
                )
            files[fname] = {
                "url": url,
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
        p = cache / (name + ".TAB" if name == "PS1P01" else (name + ".TAB").lower())
        tbl = np.loadtxt(p, delimiter=",")
        assert np.all(np.diff(tbl[:, 0]) > 0) and np.min(tbl[:, 3]) > -9
        radii = np.arange(tbl[0, 0], tbl[-1, 0] + spacing, spacing)
        # One four-bin gap in the Aldebaran profile. Interpolate in optical depth.
        tau = np.maximum(0, np.interp(radii, tbl[:, 0], tbl[:, 3]))
        profiles.append((radii[0], spacing, tau))
    states = {i: search.state_table(cache, i)[0] for i in [609, 699, 10, 399]}
    uv = disk(samples)
    outdir.mkdir(exist_ok=True, parents=True)
    with (outdir / "phoebe-lightcurve.csv").open("w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(
            [
                "earth_utc",
                "hours_from_mid",
                "planet_only",
                "planet_and_rings",
                "with_limb_darkening",
                "cassini_aldebaran",
                "cassini_gacrux",
            ]
        )
        for minutes in range(-32 * 60, 32 * 60 + 1, 5):
            t = MID + minutes / 1440
            writer.writerow(
                [
                    search.utc(search.earth_time(states, 609, t)),
                    round(minutes / 60, 4),
                    *[round(x, 6) for x in sunlight(states, t, uv, profiles)],
                ]
            )
    (outdir / "ring-sources.json").write_text(json.dumps(files, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=search.DIRECTORY)
    p.add_argument(
        "--samples",
        type=int,
        default=14,
        help="Base-two logarithm of solar disk ray count",
    )
    args = p.parse_args()
    run(args.cache, args.outdir, args.samples)
