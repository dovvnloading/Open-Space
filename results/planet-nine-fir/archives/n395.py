"""Matched-filter photometry in the 720 s, 395 nm DECam frame on the predicted tracks.

The single narrow-band visit retains its native ADU scale. A separate Legacy g template
provides a stationary-source veto, never a narrow-band flux prediction.
"""

from decam_difference import (
    HERE,
    p,
    register,
    sky_subtract,
    psf_stars,
    mad,
    sample,
    placed,
)
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as nd, signal
import numpy as np
import json, csv
from pathlib import Path


def run(src, out):
    key = "recent_n395"
    rng = np.random.default_rng(50723)
    states = p.family(
        json.loads(Path(p.__file__).with_name("orbit_input.json").read_text())
    )
    hdr = fits.getheader(src / (key + ".response"))
    a, d = p.sky(states, p.epoch(hdr["DATE-OBS"]) + hdr["EXPTIME"] / 172800)
    n = len(a)
    sig = np.full(n, np.nan)
    tests = []
    ids = rng.choice(n, 320, replace=False)
    levels = [80, 150, 250, 400, 700, 1200]
    hits = np.zeros(len(levels), int)
    for ccd in [15, 21]:
        img, h = fits.getdata(src / (key + ".response"), ccd, header=True)
        dq, hm = fits.getdata(src / (key + "_mask.response"), ccd, header=True)
        tpl, th = fits.getdata(src / f"dr11_n395ref_{ccd}.response", header=True)
        assert h["CCDNUM"] == hm["CCDNUM"] == ccd
        w = WCS(h)
        tw = WCS(th)
        x, y = w.world_to_pixel_values(a, d)
        ins = (x >= 40) & (y >= 40) & (x < img.shape[1] - 40) & (y < img.shape[0] - 40)
        xs = np.arange(
            max(0, int(x[ins].min() - 140)), min(img.shape[1], int(x[ins].max() + 140))
        )
        ys = np.arange(
            max(0, int(y[ins].min() - 140)), min(img.shape[0], int(y[ins].max() + 140))
        )
        x -= xs[0]
        y -= ys[0]
        bad = dq[np.ix_(ys, xs)] != 0
        im = sky_subtract(np.where(bad, np.nan, img[np.ix_(ys, xs)]))
        im = np.nan_to_num(im)
        tb = register(tpl, tw, w, xs, ys, 1)
        bad |= ~np.isfinite(tb)
        tb = np.nan_to_num(tb)
        sm = nd.gaussian_filter(tb, 2)
        fixed = nd.maximum_filter(sm - np.median(sm) > 5 * mad(sm), 23)
        sharp = im - nd.median_filter(im, 3) > 8 * mad(im[~bad]) + 0.35 * np.maximum(
            im, 0
        )
        mask = nd.maximum_filter(bad, 29) | nd.maximum_filter(sharp, 29) | fixed
        mask[:32] = mask[-32:] = True
        mask[:, :32] = mask[:, -32:] = True
        stars, ps = psf_stars(im, im, bad)
        kernel = ps[0]
        norm = (kernel * kernel).sum()
        fm = signal.fftconvolve(im, kernel[::-1, ::-1], mode="same") / norm
        zero = np.median(fm[~mask])
        noise = mad(fm[~mask])
        good = ins & (sample(mask.astype("f4"), x, y) == 0)
        sig[good] = (sample(fm, x[good], y[good]) - zero) / noise
        ng = 0
        for j in ids:
            if not good[j]:
                continue
            ng += 1
            ix, iy = int(round(x[j])), int(round(y[j]))
            z = im[iy - 25 : iy + 26, ix - 25 : ix + 26]
            xx = x[j] - ix + 25
            yy = y[j] - iy + 25
            for l, flux in enumerate(levels):
                stamp = np.zeros(z.shape)
                placed(stamp, (xx, yy), kernel, flux)
                injected = z + rng.poisson(np.clip(stamp * 3.8, 0, None)) / 3.8
                mf = (
                    signal.fftconvolve(injected, kernel[::-1, ::-1], mode="same") / norm
                )
                hits[l] += (sample(mf, [xx], [yy])[0] - zero) / noise >= 5.5
        tests.append(
            dict(
                ccd=ccd,
                n_psf_stars=len(stars),
                n_geom=int(ins.sum()),
                valid=int(good.sum()),
                n_plant_valid=ng,
                noise_adu=float(noise),
            )
        )
    inds = np.where(sig >= 5.5)[0]
    res = dict(
        date_obs=hdr["DATE-OBS"],
        exptime=hdr["EXPTIME"],
        filter=hdr["FILTER"],
        grid=n,
        n_valid=int(np.isfinite(sig).sum()),
        n_above=len(inds),
        max_snr=float(np.nanmax(sig)),
        ccds=tests,
        injections=[
            dict(
                flux_adu=f,
                draws=len(ids),
                valid=sum(x["n_plant_valid"] for x in tests),
                recovered=int(v),
            )
            for f, v in zip(levels, hits)
        ],
    )
    print(res)
    (out / "modern_n395_result.json").write_text(json.dumps(res, indent=2) + "\n")
    with (out / "modern_n395_flags.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["track", "ra_deg", "dec_deg", "snr"])
        wr.writerows((j, a[j], d[j], sig[j]) for j in inds)


if __name__ == "__main__":
    run(HERE / "receipts", HERE / "reduced")
