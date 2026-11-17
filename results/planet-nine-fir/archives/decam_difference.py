"""Two-night DECam image subtraction, signed moving-source filter and planted-planet recovery.

Run after obtaining each hashed input listed in decam_inputs.json. Requires NumPy,
SciPy, Astropy and the adjacent orbit propagation implementation. No precomputed
reference image enters the difference. The DR11 template checks photometry only.
"""

import argparse, csv, json
from pathlib import Path
import numpy as np
from scipy import ndimage as nd, signal
from astropy.io import fits
from astropy.wcs import WCS
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "orbits"))
import propagate as p

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(192721)
EXPOSURES = [605285, 605571]


def mad(x):
    return np.nanmedian(np.abs(x - np.nanmedian(x))) * 1.4826022185


def sky_subtract(img):
    # The two CCDs have slowly varying sky. Interpolate median values of 32px cells.
    ny, nx = img.shape
    step = 32
    by, bx = int(np.ceil(ny / step)), int(np.ceil(nx / step))
    padded = np.pad(img, ((0, by * step - ny), (0, bx * step - nx)), mode="edge")
    bg = np.nanmedian(padded.reshape(by, step, bx, step), axis=(1, 3))
    bg = np.where(np.isfinite(bg), bg, np.nanmedian(img))
    iy, ix = np.indices(img.shape)
    return img - nd.map_coordinates(
        bg, [(iy + 0.5) / step - 0.5, (ix + 0.5) / step - 0.5], order=1, mode="nearest"
    )


def register(data, src, ref, xs, ys, order):
    out = np.zeros((len(ys), len(xs)), dtype="f4")
    for i in range(0, len(ys), 80):
        xx, yy = np.meshgrid(xs, ys[i : i + 80])
        ra, dec = ref.pixel_to_world_values(xx, yy)
        x, y = src.world_to_pixel_values(ra, dec)
        out[i : i + len(xx)] = nd.map_coordinates(
            data, [y, x], order=order, mode="constant", cval=np.nan
        )
    return out


def psf_stars(a, b, invalid):
    sm = nd.gaussian_filter(b, 1.5)
    loc = (
        (sm == nd.maximum_filter(sm, size=23))
        & (sm > 45)
        & (sm < 4800)
        & (~nd.maximum_filter(invalid, size=39))
    )
    loc[:24] = loc[-24:] = False
    loc[:, :24] = loc[:, -24:] = False
    stars = []
    profiles = [[], []]
    R = 12
    yy, xx = np.mgrid[-R : R + 1, -R : R + 1]
    rad = np.hypot(xx, yy)
    for iy, ix in zip(*np.where(loc)):
        pair = []
        for sky in (a, b):
            cut = sky[iy - R : iy + R + 1, ix - R : ix + R + 1].copy()
            cut -= np.median(cut[rad > 10])
            flux = cut[rad <= 9].sum()
            if flux <= 500:
                break
            wm = np.maximum(0, cut) * (rad < 6)
            dx = (wm * xx).sum() / wm.sum()
            dy = (wm * yy).sum() / wm.sum()
            if abs(dx) > 1.5 or abs(dy) > 1.5:
                break
            centered = nd.shift(cut, [-dy, -dx], order=3, mode="constant")
            pair.append((centered, flux, dx, dy))
        if len(pair) != 2:
            continue
        # Exclude asymmetric/blended objects and galaxies with extended profiles.
        if any(
            np.sum(np.maximum(q[0], 0) * (rad < 8) * rad**2)
            / np.sum(np.maximum(q[0], 0) * (rad < 8))
            > 18
            for q in pair
        ):
            continue
        stars.append(
            [
                ix,
                iy,
                *[q[1] for q in pair],
                pair[1][2] - pair[0][2],
                pair[1][3] - pair[0][3],
            ]
        )
        for buf, q in zip(profiles, pair):
            buf.append(q[0] / q[1])
    stars = np.array(stars)
    kernels = []
    for arr in profiles:
        prof = np.median(arr, axis=0)
        prof = np.maximum(prof, 0) * (rad <= 11)
        kernels.append(prof / prof.sum())
    return stars, kernels


def sample(a, x, y):
    return nd.map_coordinates(a, [y, x], order=1, mode="constant", cval=np.nan)


def placed(inn, xy, psf, flux):
    x, y = xy
    ix, iy = int(np.rint(x)), int(np.rint(y))
    hh = psf.shape[0] // 2
    rr = nd.shift(psf, (y - iy, x - ix), order=3)
    rr /= rr.sum()
    inn[iy - hh : iy + hh + 1, ix - hh : ix + hh + 1] += flux * rr


def lobes(fa, fb, noise):
    # A constant-brightness mover leaves opposite-sign PSF residuals at its two
    # predicted positions. Reject one-sided features and >4:1 lobe ratios.
    return (np.minimum(fa, -fb) >= 1.75 * noise / np.sqrt(2)) & (
        np.minimum(fa, -fb) >= 0.25 * np.maximum(fa, -fb)
    )


def filtered(a, b, psfs, sc):
    pa, pb = psfs
    d = signal.fftconvolve(a, pb, mode="same") - sc * signal.fftconvolve(
        b, pa, mode="same"
    )
    pk = signal.fftconvolve(pa, pb, mode="full")
    return signal.fftconvolve(d, pk[::-1, ::-1], mode="same") / np.sum(pk * pk)


def validate_dense_bank(
    a,
    b,
    psfs,
    scale,
    wcs0,
    xs,
    ys,
    states,
    times,
    bad,
    scoremap,
    region,
    vecs,
    noise,
    center,
    threshold,
    zp,
):
    # Independent discretization, randomly perturbed in sky position by a uniform-area
    # 12-arcsec disk. This is a sensitivity experiment, not an astrometric prior.
    rng = np.random.default_rng(579921)
    teststates = p.family(
        json.loads(Path(p.__file__).with_name("orbit_input.json").read_text()),
        13,
        32,
        56,
    )
    inds = rng.choice(len(teststates[0]), size=400, replace=False)
    off = rng.uniform(0, 2 * np.pi, len(inds))
    rho = 12 * np.sqrt(rng.uniform(size=len(inds)))
    positions = []
    for t in times:
        ra, dec = p.sky(teststates, t)
        ra = ra[inds] + rho * np.cos(off) / 3600 / np.cos(np.deg2rad(dec[inds]))
        dec = dec[inds] + rho * np.sin(off) / 3600
        x, y = wcs0.world_to_pixel_values(ra, dec)
        positions.append((x - xs[0], y - ys[0]))
    mags = [21.0, 22.0, 22.5, 23.0, 23.25, 23.5, 23.75, 24.0]
    counts = np.zeros(len(mags), int)
    tested = 0
    local = np.mgrid[-2:3, -2:3]
    ox = local[1].ravel()
    oy = local[0].ravel()
    for j in range(len(inds)):
        x = np.array([v[0][j] for v in positions])
        y = np.array([v[1][j] for v in positions])
        ix, iy = np.rint([x[0], y[0]]).astype(int)
        if not region[iy, ix]:
            continue
        xbase = ix + ox
        ybase = iy + oy
        close = (xbase - x[0]) ** 2 + (ybase - y[0]) ** 2 <= 2.5**2
        xbase = xbase[close]
        ybase = ybase[close]
        xb = xbase[:, None] + vecs[:, 0]
        yb = ybase[:, None] + vecs[:, 1]
        valid = (~bad[ybase, xbase])[:, None] & (~bad[yb, xb])
        if not np.any(valid):
            continue
        tested += 1
        bx = ix - 50
        by = iy - 50
        x = x - bx
        y = y - by
        sl = np.s_[by : by + 101, bx : bx + 101]
        rawa = a[sl]
        rawb = b[sl]
        for mi, mag in enumerate(mags):
            aa = rawa.copy()
            bb = rawb.copy()
            flux = 10 ** (0.4 * (zp - mag))
            for k, (im, star, nadu) in enumerate(
                zip((aa, bb), psfs, (flux, flux / scale))
            ):
                model = np.zeros(im.shape)
                placed(model, (x[k], y[k]), star, nadu)
                im += rng.poisson(np.clip(model * 3.83, 0, None)) / 3.83
            fm = filtered(aa, bb, psfs, scale)
            sv = (
                fm[ybase - by, xbase - bx][:, None] - fm[yb - by, xb - bx] - center
            ) / noise
            counts[mi] += np.any(
                (sv >= threshold)
                & valid
                & lobes(
                    fm[ybase - by, xbase - bx][:, None], fm[yb - by, xb - bx], noise
                )
            )
    return dict(
        n_draws=len(inds),
        n_valid=tested,
        perturbation_radius_arcsec=12,
        grid=[13, 32, 56],
        localization_radius_pix=2.5,
        source_poisson_gain_e_per_adu=3.83,
        completeness=[
            dict(
                r_ab=m,
                n_recovered=int(v),
                fraction_all=float(v / len(inds)),
                fraction_eligible=float(v / tested),
            )
            for m, v in zip(mags, counts)
        ],
    )


def main(inputs, outdir):
    outdir.mkdir(exist_ok=True, parents=True)
    states = p.family(
        json.loads(Path(p.__file__).with_name("orbit_input.json").read_text())
    )
    coords = []
    data = []
    masks = []
    headers = []
    wcs = []
    stamps = []
    for ident in EXPOSURES:
        img, hdr = fits.getdata(inputs / f"DECam_{ident}.response", header=True)
        mask, mhdr = fits.getdata(inputs / f"mask_{ident}.response", header=True)
        assert (
            hdr["EXPNUM"] == ident == mhdr["EXPNUM"] and hdr["CCDNUM"] == mhdr["CCDNUM"]
        )
        w = WCS(hdr)
        t = p.epoch(hdr["DATE-OBS"]) + hdr["EXPTIME"] / 172800
        coords.append(p.sky(states, t))
        wcs.append(w)
        headers.append(hdr)
        stamps.append(t)
        data.append(img.astype("f4"))
        masks.append(mask != 0)
    xy = [np.array(wcs[0].world_to_pixel_values(*c)) for c in coords]
    box = np.hstack(xy)
    xs = np.arange(int(np.floor(box[0].min())) - 370, int(np.ceil(box[0].max())) + 371)
    ys = np.arange(int(np.floor(box[1].min())) - 370, int(np.ceil(box[1].max())) + 371)
    xs = xs[(xs >= 0) & (xs < data[0].shape[1])]
    ys = ys[(ys >= 0) & (ys < data[0].shape[0])]
    xy = [v - np.array([[xs[0]], [ys[0]]]) for v in xy]
    a = data[0][np.ix_(ys, xs)]
    b = register(data[1], wcs[1], wcs[0], xs, ys, 1)
    m = (
        masks[0][np.ix_(ys, xs)]
        | (register(masks[1].astype("f4"), wcs[1], wcs[0], xs, ys, 1) != 0)
        | ~np.isfinite(a + b)
    )
    a = sky_subtract(np.where(m, np.nan, a))
    b = sky_subtract(np.where(m, np.nan, b))
    a = np.nan_to_num(a)
    b = np.nan_to_num(b)
    stars, psfs = psf_stars(a, b, m)
    dx, dy = np.median(stars[:, 4:], axis=0)
    # Refine relative translation using the stationary stars. Retain A's absolute WCS.
    b = nd.shift(b, [-dy, -dx], order=1, mode="nearest")
    m |= nd.shift(m.astype("f4"), [-dy, -dx], order=1, mode="nearest") > 0
    stars, psfs = psf_stars(a, b, m)
    scale = np.median(stars[:, 2] / stars[:, 3])
    starrat = mad(np.log10(stars[:, 2] / stars[:, 3])) * 2.5
    template, hh = fits.getdata(inputs / "template_605571.response", header=True)
    tpl = register(template, WCS(hh), wcs[0], xs, ys, 1)
    zps = []
    radius = 10
    u, v = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    rr = np.hypot(u, v)
    for j in range(2):
        vals = []
        for sx, sy, fa, fb, *_ in stars:
            x, y = int(sx), int(sy)
            stamp = tpl[y - radius : y + radius + 1, x - radius : x + radius + 1]
            ft = (stamp - np.nanmedian(stamp[rr > 9]))[rr < 8].sum()
            f = [fa, fb][j]
            if ft > 0 and np.isfinite(ft):
                vals.append(22.5 + 2.5 * np.log10(f / ft))
        zps.append((float(np.median(vals)), float(mad(vals)), len(vals)))
    # Reprojection preserves the near-identical native solid angles of consecutive S12 exposures.
    # The archive's generic Gaia MAGZERO differs from the field-based Legacy r calibration.
    zp = zps[0][0]
    scoremap = filtered(a, b, psfs, scale)
    dilated = nd.maximum_filter(m, size=29)
    # Flag neighborhoods of strong stationary sources: interpolation/subtraction residual floor.
    bright = nd.maximum_filter(np.minimum(a, b * scale), size=27) > 220

    # Unflagged particle tracks in the first exposure created narrow, curved features.
    # A resolved stellar PSF is much broader than these single-pixel gradients.
    def sharp(im):
        med = nd.median_filter(im, size=3)
        return (im - med) > 8 * mad(im[~m]) + 0.35 * np.maximum(im, 0)

    particle = nd.maximum_filter(sharp(a) | sharp(b), size=37)
    bad = dilated | bright | particle
    x0, y0 = xy[0]
    x1, y1 = xy[1]
    ok = (sample(bad.astype("f4"), x0, y0) == 0) & (
        sample(bad.astype("f4"), x1, y1) == 0
    )
    raw = sample(scoremap, x0, y0) - sample(scoremap, x1, y1)
    n = len(raw)
    # Control dipoles at displaced positions, same displacement and detector neighborhoods.
    jit = RNG.uniform(-100, 100, (2, 32000))
    sel = RNG.integers(n, size=jit.shape[1])
    ex, ey = jit
    cx0, cy0 = x0[sel] + ex, y0[sel] + ey
    cx1, cy1 = x1[sel] + ex, y1[sel] + ey
    cgood = (sample(bad.astype("f4"), cx0, cy0) == 0) & (
        sample(bad.astype("f4"), cx1, cy1) == 0
    )
    cn = (sample(scoremap, cx0, cy0) - sample(scoremap, cx1, cy1))[cgood]
    noise = mad(cn)
    center = np.median(cn)
    score = (raw - center) / noise
    threshold = max(5.5, float(np.quantile((cn - center) / noise, 0.9995)))
    survive = (
        ok
        & (score >= threshold)
        & lobes(sample(scoremap, x0, y0), sample(scoremap, x1, y1), noise)
    )
    # Dense trajectory-agnostic check: evaluate every pixel in the box around both epoch predictions,
    # and every rounded relative motion vector seen in the sampled grid. Expand the bounds by
    # 290 pixels (~76 arcsec), clipped to the usable CCD interior. Unit-pixel positioning has maximum quantization of 0.19 arcsec per epoch.
    vecs = np.unique(np.rint(np.stack((x1 - x0, y1 - y0), axis=-1)).astype(int), axis=0)
    yy, xx = np.indices(scoremap.shape)
    margin = 290
    region = (
        (xx >= 35)
        & (xx < xx.shape[1] - 35)
        & (yy >= 35)
        & (yy < xx.shape[0] - 35)
        & (xx >= x0.min() - margin)
        & (xx <= x0.max() + margin)
        & (yy >= y0.min() - margin)
        & (yy <= y0.max() + margin)
    )
    eligible = np.zeros(scoremap.shape, dtype=bool)
    best = np.full(scoremap.shape, -np.inf)
    vbest = np.zeros(scoremap.shape, int)
    for i, (vx, vy) in enumerate(vecs):
        ds = (scoremap - np.roll(scoremap, (-vy, -vx), axis=(0, 1)) - center) / noise
        valid = region & ~bad & ~np.roll(bad, (-vy, -vx), axis=(0, 1))
        eligible |= valid
        fs = np.roll(scoremap, (-vy, -vx), axis=(0, 1))
        improves = valid & lobes(scoremap, fs, noise) & (ds > best)
        best[improves] = ds[improves]
        vbest[improves] = i
    peaks = (best == nd.maximum_filter(best, size=13)) & (best >= threshold)
    py, px = np.where(peaks)
    rr, dd = wcs[0].pixel_to_world_values(xs[0] + px, ys[0] + py)
    order = np.argsort(best[peaks])[::-1]
    nearby = []
    for jj in order:
        k = np.argmin((x0 - px[jj]) ** 2 + (y0 - py[jj]) ** 2)
        nearby.append(
            dict(
                ra=float(rr[jj]),
                dec=float(dd[jj]),
                score=float(best[py[jj], px[jj]]),
                dx=int(vecs[vbest[py[jj], px[jj]], 0]),
                dy=int(vecs[vbest[py[jj], px[jj]], 1]),
                nearest_track=int(k),
                separation_pix=float(np.hypot(x0[k] - px[jj], y0[k] - py[jj])),
            )
        )
    (outdir / "dense_peaks.json").write_text(json.dumps(nearby, indent=2) + "\n")
    print("dense candidates", len(nearby), nearby[:3], flush=True)
    offgrid = validate_dense_bank(
        a,
        b,
        psfs,
        scale,
        wcs[0],
        xs,
        ys,
        states,
        stamps,
        bad,
        scoremap,
        region,
        vecs,
        noise,
        center,
        threshold,
        zp,
    )
    (outdir / "dense_injections.json").write_text(json.dumps(offgrid, indent=2) + "\n")
    from planet_injections import measure

    planets = measure(
        p,
        placed,
        filtered,
        lobes,
        a,
        b,
        psfs,
        scale,
        wcs[0],
        xs,
        ys,
        stamps,
        bad,
        region,
        vecs,
        noise,
        center,
        threshold,
        zp,
    )
    (outdir / "model_recovery.json").write_text(json.dumps(planets, indent=2) + "\n")
    # Plant one model at a time into copies of the registered science pixels, before PSF
    # matching/subtraction/filtering. Every position is drawn uniformly from finite grid rows.
    # Epoch displacements are preserved. Avoid cross-talk between simultaneous injections.
    mags = np.arange(21.0, 25.01, 0.5)
    plant_ids = RNG.choice(n, size=320, replace=False)
    records = []
    by_mag = []
    for mag in mags:
        recovered = 0
        tested = 0
        for j in plant_ids:
            if not ok[j]:
                continue
            bx = int(round(0.5 * (x0[j] + x1[j]))) - 48
            by = int(round(0.5 * (y0[j] + y1[j]))) - 48
            sl = np.s_[by : by + 97, bx : bx + 97]
            x = np.array([x0[j], x1[j]]) - bx
            y = np.array([y0[j], y1[j]]) - by
            aa = a[sl].copy()
            bb = b[sl].copy()
            flux = 10 ** (0.4 * (zp - mag))
            placed(aa, (x[0], y[0]), psfs[0], flux)
            placed(bb, (x[1], y[1]), psfs[1], flux / scale)
            fm = filtered(aa, bb, psfs, scale)
            det = (sample(fm, x, y) @ [1.0, -1.0] - center) / noise
            tested += 1
            hit = det >= threshold and lobes(*sample(fm, x, y), noise)
            recovered += int(hit)
            records.append([j, mag, float(det), int(hit)])
        by_mag.append(
            dict(
                r_ab=float(mag),
                planted_unmasked=tested,
                recovered=recovered,
                grid_draws=len(plant_ids),
                fraction_all_draws=recovered / len(plant_ids),
                fraction_valid=recovered / tested,
            )
        )
    result = dict(
        dense_max_snr=float(best[np.isfinite(best)].max()),
        dense_area_pixels=int(region.sum()),
        dense_eligible_pixels=int(eligible.sum()),
        dense_trials_motion_vectors=vecs.tolist(),
        dense_peaks=len(nearby),
        exposure_ids=EXPOSURES,
        mjd_mid=[float(h["MJD-OBS"] + h["EXPTIME"] / 172800) for h in headers],
        orbit_grid=n,
        mask_valid=int(ok.sum()),
        max_snr=float(score[ok].max()),
        n_above=int(survive.sum()),
        threshold=threshold,
        empirical_control_noise_adu=float(noise),
        controls=len(cn),
        control_median_adu=float(center),
        empirical_shift_xy_pix=[float(dx), float(dy)],
        n_psf_stars=len(stars),
        scale_b_to_a=float(scale),
        relative_scatter_mag=float(starrat),
        local_legacy_r_zp_adu=zps,
        peak_track_index=int(np.where(ok, score, -np.inf).argmax()),
        median_track_displacement_pix=float(np.median(np.hypot(x1 - x0, y1 - y0))),
        injections=by_mag,
    )
    (outdir / "difference_result.json").write_text(json.dumps(result, indent=2) + "\n")
    for name, fields, rows in [
        (
            "difference_tracks",
            ["track", "ra_a", "dec_a", "ra_b", "dec_b", "valid", "snr", "flagged"],
            zip(
                range(n),
                *coords[0],
                *coords[1],
                ok.astype(int),
                score,
                survive.astype(int),
            ),
        ),
        ("difference_injections", ["track", "r_AB", "snr", "recovered"], records),
    ]:
        import gzip

        with gzip.open(outdir / (name + ".csv.gz"), "wt", newline="") as f:
            w = csv.writer(f)
            w.writerow(fields)
            w.writerows(rows)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", type=Path, default=HERE / "receipts")
    ap.add_argument("--out", type=Path, default=HERE / "reduced")
    args = ap.parse_args()
    main(args.inputs, args.out)
