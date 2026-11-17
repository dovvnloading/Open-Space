"""Pixel injection of reflected-light models at their distances and propagated coordinates."""

import json
from pathlib import Path
import numpy as np
from scipy import ndimage as nd
from scipy.stats import binomtest


def distance_at(p, s, t):
    # Observer at reception, object and Sun at emission.
    r, v, tref, _ = s
    obs = p.earth(t)
    light = p.norm(r - obs) / p.C
    for _ in range(3):
        loc = p.kepler(r, v, t - tref - light)
        delta = p.norm(loc - obs)
        light = delta / p.C
    hel, bary = p.observers(t - light)
    return p.norm(loc - (bary - hel)), delta


def offset(coord, east, north, p):
    n = p.unit(*coord)
    ra = np.deg2rad(coord[0])
    e = np.array([-np.sin(ra), np.cos(ra), 0])
    q = n + (east * e + north * np.cross(n, e)) / p.ARCSEC
    q /= p.norm(q)
    return [np.rad2deg(np.arctan2(q[1], q[0])) % 360, np.rad2deg(np.arcsin(q[2]))]


def population(p, wcs0, xs, ys, times, stress, rng, n=400):
    inputs = json.loads(Path(p.__file__).with_name("orbit_input.json").read_text())
    pa = np.deg2rad(79)
    iras = [
        (0, 0),
        (18 * np.sin(pa), 18 * np.cos(pa)),
        (-18 * np.sin(pa), -18 * np.cos(pa)),
        (7 * np.cos(pa), -7 * np.sin(pa)),
        (-7 * np.cos(pa), 7 * np.sin(pa)),
    ]
    akari = [(0, 0), (40, 0), (-40, 0), (0, 40), (0, -40)]
    pairs = [(i, j) for i in iras for j in akari] if stress else [((0, 0), (0, 0))]
    choices = rng.integers(len(pairs), size=n)
    out = []
    for k, (ii, aa) in enumerate(pairs):
        num = np.count_nonzero(choices == k)
        if not num:
            continue
        inp = dict(
            inputs,
            iras_radec_deg=offset(inputs["iras_radec_deg"], *ii, p),
            akari_radec_deg=offset(inputs["akari_radec_deg"], *aa, p),
        )
        s = p.family(inp, 13, 32, 56)
        j = rng.integers(len(s[0]), size=num)
        ss = tuple([x[j] for x in s[:3]] + [s[3]])
        rd = distance_at(p, ss, times[0])
        pos = []
        for t in times:
            a, d = p.sky(ss, t)
            x, y = wcs0.world_to_pixel_values(a, d)
            pos.extend([x - xs[0], y - ys[0]])
        for item in zip(*pos, *rd, *ss[0].T, *ss[1].T, ss[2]):
            out.append(item)
    return np.array(out)


def measure(
    p,
    placed,
    filtered,
    lobes,
    a,
    b,
    psfs,
    scale,
    wcs0,
    xs,
    ys,
    times,
    bad,
    region,
    vecs,
    noise,
    center,
    threshold,
    zp,
    save_draws_path=None,
):
    rng = np.random.default_rng(813481)
    models = [
        (mass, radius, albedo)
        for mass, radius in [
            (5, 2.71),
            (5, 2.94),
            (5, 4.12),
            (10, 3.16),
            (10, 3.46),
            (10, 3.70),
            (10, 5.09),
        ]
        for albedo in [0.10, 0.30, 0.75]
    ]
    outputs = []
    draw_rows = []
    for stress in [False, True]:
        rows = population(p, wcs0, xs, ys, times, stress, rng)
        rr, dd = rows[:, 4:6].T
        hits = np.zeros(len(models), int)
        trial_valid = np.zeros(len(models), int)
        mags = []
        eligible = 0
        recovered_each = np.zeros((len(rows), len(models)), bool)
        for m, R, A in models:
            mags.append(
                -26.93
                - 2.5 * np.log10(A * (R * 6371000.0 / 149597870700.0 / (rr * dd)) ** 2)
            )
        for j, (x0, y0, x1, y1, r, delta, *state) in enumerate(rows):
            ix, iy = np.rint([x0, y0]).astype(int)
            bx, by = ix - 48, iy - 48
            if (
                bx < 0
                or by < 0
                or by + 97 > a.shape[0]
                or bx + 97 > a.shape[1]
                or not region[iy, ix]
            ):
                continue
            oy, ox = np.mgrid[-2:3, -2:3]
            qx, qy = ix + ox.ravel(), iy + oy.ravel()
            near = np.hypot(qx - x0, qy - y0) <= 2.5
            qx = qx[near]
            qy = qy[near]
            nx = qx[:, None] + vecs[:, 0]
            ny = qy[:, None] + vecs[:, 1]
            valid = (~bad[qy, qx])[:, None] & (~bad[ny, nx])
            if not valid.any():
                continue
            eligible += 1
            sl = np.s_[by : by + 97, bx : bx + 97]
            ra = a[sl]
            rb = b[sl]
            for mi, ms in enumerate(mags):
                flux = 10 ** (0.4 * (zp - ms[j]))
                ii = ra.copy()
                jj = rb.copy()
                for im, psf, x, y, F in zip(
                    [ii, jj], psfs, [x0, x1], [y0, y1], [flux, flux / scale]
                ):
                    star = np.zeros(im.shape)
                    placed(star, (x - bx, y - by), psf, F)
                    im += rng.poisson(np.clip(star * 3.83, 0, None)) / 3.83
                # Reapply data-dependent bright and particle masks to the planted stamps.
                stat = nd.maximum_filter(np.minimum(ii, jj * scale), size=27) > 220
                sh = (
                    ii - nd.median_filter(ii, 3)
                    > 8 * np.median(np.abs(ra - np.median(ra))) * 1.4826
                    + 0.35 * np.maximum(ii, 0)
                ) | (
                    jj - nd.median_filter(jj, 3)
                    > 8 * np.median(np.abs(rb - np.median(rb))) * 1.4826
                    + 0.35 * np.maximum(jj, 0)
                )
                ignore = stat | nd.maximum_filter(sh, size=37)
                ok = (
                    valid
                    & (~ignore[qy - by, qx - bx])[:, None]
                    & (~ignore[ny - by, nx - bx])
                )
                if not ok.any():
                    continue
                trial_valid[mi] += 1
                f = filtered(ii, jj, psfs, scale)
                fa = f[qy - by, qx - bx][:, None]
                fb = f[ny - by, nx - bx]
                hit = np.any(
                    ok
                    & lobes(fa, fb, noise)
                    & ((fa - fb - center) >= threshold * noise)
                )
                hits[mi] += hit
                recovered_each[j, mi] = hit
        for j, row in enumerate(rows):
            draw_rows.append(
                [
                    "centroid_stress" if stress else "published_centroids",
                    j,
                    *row[6:],
                    row[4],
                    row[5],
                    *[int(recovered_each[j, k]) for k in [0, 4, 16]],
                ]
            )
        for (m, R, A), v, tr, mag in zip(models, hits, trial_valid, mags):
            ci = binomtest(int(v), len(rows)).proportion_ci(0.95, method="wilson")
            outputs.append(
                dict(
                    experiment="centroid_stress" if stress else "published_centroids",
                    mass_Mearth=m,
                    radius_Rearth=R,
                    geometric_albedo_r=A,
                    n_drawn=len(rows),
                    n_eligible_before_plant=eligible,
                    n_eligible_after_plant=int(tr),
                    n_recovered=int(v),
                    fraction=float(v / len(rows)),
                    ci95_wilson=[ci.low, ci.high],
                    r_AB_range=[float(mag.min()), float(mag.max())],
                )
            )
    if save_draws_path is not None:
        import csv, gzip

        with gzip.open(save_draws_path, "wt", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "experiment",
                    "draw",
                    "x_au",
                    "y_au",
                    "z_au",
                    "vx_au_d",
                    "vy_au_d",
                    "vz_au_d",
                    "reference_emission_tt_days_since_j2000",
                    "des_heliocentric_au",
                    "des_geocentric_au",
                    "hit_mass5_R2.71_p0.10",
                    "hit_mass5_R2.94_p0.30",
                    "hit_mass10_R3.70_p0.30",
                ]
            )
            w.writerows(draw_rows)
    return outputs
