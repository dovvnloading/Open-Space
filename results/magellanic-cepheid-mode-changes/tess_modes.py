"""Test localization and relative overtone strengths in all seven TESS visits."""

import json
import numpy as np
from scipy.linalg import eigh
from astropy.io import fits
from astropy.wcs import WCS
from screen import ROOT, INPUT, catalog
from tess_fetch import RA, DEC

PERIODS = catalog()["OGLE-SMC-CEP-2332"][1:]


def reduce(s):
    path = INPUT / "tess" / f"s{s}.fits"
    with fits.open(path, memmap=False) as h:
        data = h[1].data
        img = np.array(data["FLUX"])
        t = np.array(data["TIME"])
        use = (
            (data["QUALITY"] == 0)
            & np.isfinite(t)
            & np.all(np.isfinite(img), axis=(1, 2))
        )
        img = img[use]
        t = t[use]
        x, y = WCS(h[2].header, naxis=2).world_to_pixel_values(RA, DEC)
        dates = [h[0].header["DATE-OBS"], h[0].header["DATE-END"]]
    yy, xx = np.indices(img.shape[1:])
    dist = np.hypot(xx - x, yy - y).ravel()
    Y = img.reshape(len(t), -1).astype(float)
    Y -= np.median(Y, axis=0)
    # Bin on 0.025-day cells to compare short/long cadences without pseudo-replication.
    key = np.floor(t / 0.025).astype(int)
    n = np.unique(key, return_inverse=True)[1]
    counts = np.bincount(n)
    mat = np.zeros((len(counts), Y.shape[1]))
    np.add.at(mat, n, Y)
    Y = mat / counts[:, None]
    t = np.bincount(n, weights=t) / counts
    ref = Y[:, dist > 3.1]
    scale = np.maximum(
        1.4826 * np.median(abs(ref - np.median(ref, axis=0)), axis=0), 0.5
    )
    ref = np.clip(ref / scale, -7, 7)
    _, v = eigh(ref.T @ ref, subset_by_index=(ref.shape[1] - 6, ref.shape[1] - 1))
    pc = ref @ v[:, ::-1]
    tc = t - np.median(t)

    def run(npc):
        terms = [np.ones(len(t)), tc / 20, (tc / 20) ** 2, *pc[:, :npc].T]
        freqs = [
            1 / PERIODS[0],
            1 / PERIODS[1],
            2 / PERIODS[0],
            2 / PERIODS[1],
            sum(1 / z for z in PERIODS),
            1 / PERIODS[1] - 1 / PERIODS[0],
        ]
        k = len(terms)
        for freq in freqs:
            phi = 2 * np.pi * tc * freq
            terms.extend([np.sin(phi), np.cos(phi)])
        X = np.column_stack(terms)
        coeff = np.linalg.lstsq(X, Y, rcond=None)[0]
        res = Y - X @ coeff
        # An aperture fixed on the OGLE position and target-centered alternative aperture.
        out = {}
        radii = [1.1, 1.6, 2.1]
        for radius in radii:
            ap = dist < radius
            # Fit after clipping isolated instrumental excursions in the aperture.
            residual = res[:, ap].sum(axis=1)
            med = np.median(residual)
            sig = 1.4826 * np.median(abs(residual - med))
            keep = abs(residual - med) < 5 * sig
            b = np.linalg.lstsq(X[keep], Y[keep][:, ap].sum(axis=1), rcond=None)[0]
            a1, a2 = (float(np.hypot(*b[j : j + 2])) for j in [k, k + 2])
            # Nearby frequency controls stay clear of both overtone frequencies and harmonics.
            yres = Y[keep][:, ap].sum(axis=1) - X[keep] @ b
            control = []
            for f in np.linspace(0.4, 3.3, 75):
                if min(abs(f - z) for z in freqs) < 0.12:
                    continue
                phi = 2 * np.pi * tc[keep] * f
                u = np.linalg.lstsq(
                    np.column_stack([np.sin(phi), np.cos(phi)]), yres, rcond=None
                )[0]
                control.append(np.hypot(*u))
            out[str(radius)] = dict(
                a1=a1,
                a2=a2,
                ratio=a2 / a1,
                control95=float(np.percentile(control, 95)),
                kept=int(sum(keep)),
            )
        # Offset of strongest amplitude pixel strictly in central 5x5; compare typical offsets to 21" pixel.
        central = dist < 2.75
        peak = []
        for j in [k, k + 2]:
            a = np.hypot(coeff[j], coeff[j + 1])
            pix = np.argmax(np.where(central, a, -1))
            peak.append(
                dict(
                    dx=float(xx.ravel()[pix] - x),
                    dy=float(yy.ravel()[pix] - y),
                    amplitude=float(a[pix]),
                    center=float(a[np.argmin(dist)]),
                )
            )
        return dict(apertures=out, peaks=peak)

    return dict(
        sector=s,
        dates=dates,
        n_cadences=int(sum(use)),
        n_bins=len(t),
        tHJD_approx=float(t.mean() + 7000),
        x=float(x),
        y=float(y),
        models={str(i): run(i) for i in [2, 4, 6]},
    )


if __name__ == "__main__":
    results = [reduce(s) for s in [1, 27, 28, 67, 68, 94, 95]]
    (ROOT / "tess_checks.json").write_text(json.dumps(results, indent=2) + "\n")
    for r in results:
        v = r["models"]["4"]
        print(
            r["sector"],
            r["dates"][0][:10],
            [
                (round(z, 2))
                for z in [
                    v["apertures"]["1.6"][k] for k in ["a1", "a2", "ratio", "control95"]
                ]
            ],
            [(round(p["dx"], 1), round(p["dy"], 1)) for p in v["peaks"]],
        )
