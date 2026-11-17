"""Two November 2024 g exposures versus a pre-September-2024 Legacy DR11 g template.

Template-subtracted photometry tests the model directions on one night. The 62 s
baseline does not measure the expected motion. CCD 50 alone is used for photometry.
"""

from decam_difference import (
    HERE,
    p,
    sky_subtract,
    register,
    psf_stars,
    filtered,
    sample,
    placed,
    mad,
)
import csv, json
from pathlib import Path
from datetime import timedelta
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy import ndimage as nd


def measure(src, dest):
    rng = np.random.default_rng(195731)
    states = p.family(
        json.loads(Path(p.__file__).with_name("orbit_input.json").read_text())
    )
    keys = ["recent_g1", "recent_g2"]
    ccd = 50
    template, th = fits.getdata(src / "dr11_g_50.response", header=True)
    tw = WCS(th)
    h1 = fits.getheader(src / "recent_g1.response", ccd)
    ref = WCS(h1)
    model = []
    results = []
    expmap = []
    cutm = []
    cal = []
    allcoords = []
    filters = []
    for key in keys:
        h = fits.getheader(src / (key + ".response"))
        coords = p.sky(states, p.epoch(h["DATE-OBS"]) + h["EXPTIME"] / 172800)
        allcoords.append(coords)
        model.append(np.array(ref.world_to_pixel_values(*coords)))
    xy = np.hstack(model)
    xs = np.arange(
        max(15, int(xy[0].min() - 150)), min(h1["NAXIS1"] - 15, int(xy[0].max() + 150))
    )
    ys = np.arange(
        max(15, int(xy[1].min() - 150)), min(h1["NAXIS2"] - 15, int(xy[1].max() + 150))
    )
    model = [q - np.array([[xs[0]], [ys[0]]]) for q in model]
    tb = register(template, tw, ref, xs, ys, 1)
    tm = ~np.isfinite(tb)
    tb = np.nan_to_num(sky_subtract(tb))
    smooth = nd.gaussian_filter(tb, 2)
    static = nd.maximum_filter(smooth > 5 * mad(smooth), size=23)
    for k, key in enumerate(keys):
        im, h = fits.getdata(src / (key + ".response"), ccd, header=True)
        m, hh = fits.getdata(src / (key + "_mask.response"), ccd, header=True)
        assert h["CCDNUM"] == hh["CCDNUM"] == ccd
        w = WCS(h)
        a = register(im, w, ref, xs, ys, 1)
        bad = tm | (register(m.astype("f4"), w, ref, xs, ys, 1) != 0) | ~np.isfinite(a)
        a = np.nan_to_num(sky_subtract(np.where(bad, np.nan, a)))
        s0 = 10 ** (0.4 * (fits.getheader(src / (key + ".response"))["MAGZERO"] - 22.5))
        b = tb * s0
        stars, ps = psf_stars(a, b, bad)
        dx, dy = np.median(stars[:, 4:], axis=0)
        b = nd.shift(b, [-dy, -dx], order=1)
        stars, ps = psf_stars(a, b, bad)
        scale = np.median(stars[:, 2] / stars[:, 3])
        zp = 22.5 + 2.5 * np.log10(s0 * scale)
        f = filtered(a, b, ps, scale)
        # in exposure ADU
        msk = static | nd.maximum_filter(bad, size=29)
        msk[:32] = msk[-32:] = True
        msk[:, :32] = msk[:, -32:] = True
        e = mad(f[~msk])
        fmed = np.median(f[~msk])
        fm = (f - fmed) / e
        results.append(
            dict(
                exposure=key,
                n_psf_stars=len(stars),
                sigma_adu=e,
                zp_ab=zp,
                relative_scale=scale,
                fit_offsets_pix=[dx, dy],
            )
        )
        expmap.append(fm)
        cutm.append(msk)
        cal.append((a, b, ps, scale, e, fmed, zp))
    msk = cutm[0] | cutm[1]
    avg = np.mean(expmap, axis=0)
    jn = mad(avg[~msk])
    scores = np.stack([sample(fm, *q) for fm, q in zip(expmap, model)])
    valid = np.logical_and.reduce(
        [sample(m.astype("f4"), *q) == 0 for m, q in zip(cutm, model)]
    )
    sn = np.mean(scores, axis=0) / jn
    flags = valid & (sn >= 5.5) & (scores.min(0) >= 2)
    with (dest / "modern_g_tracks.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["model", "ra", "dec", "joint_snr", "g1_snr", "g2_snr"])
        wr.writerows(
            (j, allcoords[0][0][j], allcoords[0][1][j], sn[j], *scores[:, j])
            for j in np.flatnonzero(flags)
        )
    tests = []
    ids = rng.choice(len(valid), size=320, replace=False)
    mags = [19.0, 20.0, 21.0, 21.5, 22.0, 22.5, 23.0]
    for mag in mags:
        hits = 0
        tested = 0
        for j in ids:
            if not valid[j]:
                continue
            sig = []
            for (a, b, ps, scale, e, base, zp), q in zip(cal, model):
                x, y = q[:, j]
                ix, iy = int(np.rint(x)), int(np.rint(y))
                sl = np.s_[iy - 36 : iy + 37, ix - 36 : ix + 37]
                aa = a[sl].copy()
                bb = b[sl]
                xx = x - ix + 36
                yy = y - iy + 36
                pp = np.zeros(aa.shape)
                placed(pp, (xx, yy), ps[0], 10 ** (0.4 * (zp - mag)))
                aa += rng.poisson(np.clip(pp * 3.8, 0, None)) / 3.8
                ff = filtered(aa, bb, ps, scale)
                sig.append((sample(ff, [xx], [yy])[0] - base) / e)
            tested += 1
            hits += np.mean(sig) / jn >= 5.5 and min(sig) >= 2
        tests.append(dict(g_AB=mag, draws=len(ids), valid=tested, recovered=int(hits)))
    res = dict(
        exposures=results,
        ccd=ccd,
        grid_nodes=len(valid),
        n_valid=int(valid.sum()),
        n_flagged=int(flags.sum()),
        max_snr=float(sn[valid].max()),
        joint_noise=jn,
        injections=tests,
    )
    (dest / "modern_g_result.json").write_text(json.dumps(res, indent=2) + "\n")
    print(res, flush=True)


if __name__ == "__main__":
    measure(HERE / "receipts", HERE / "reduced")
