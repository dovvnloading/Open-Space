"""Fit the known 3.599551-day Cepheid and its first four harmonics in OGLE and TESS."""

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from open_space.ffi_archive import acquire
from open_space.tess_archive import atomic_bytes, download

PERIOD = 3.599551
CENTER = (79.02466, -69.34685)
STAR = (79.03048646142, -69.34643753244)


def fit(t, flux, error, extra=None):
    angle = 2 * np.pi * (t - 7000) / PERIOD
    basis = np.column_stack(
        [np.ones(len(t))]
        + [fn(k * angle) for k in range(1, 5) for fn in (np.cos, np.sin)]
    )
    if extra is not None:
        basis = np.column_stack([basis, extra])
    beta = np.linalg.lstsq(basis / error[:, None], flux / error, rcond=None)[0]
    coeff = beta[1:9].reshape(4, 2)
    amp = np.hypot(*coeff.T)
    phase = np.arctan2(coeff[:, 1], coeff[:, 0]) / (2 * np.pi)
    return {
        "epochs": len(t),
        "semiamplitudes": amp.tolist(),
        "phase_cycles": phase.tolist(),
        "third_over_fundamental": amp[2] / amp[0],
        "third_minus_three_fundamental_phase_cycles": (phase[2] - 3 * phase[0] + 0.5)
        % 1
        - 0.5,
    }


def tess(cube):
    t, img, err, wcs = cube
    keep = np.all(np.isfinite(img) & np.isfinite(err) & (err > 0), axis=(1, 2))
    t, img, err = t[keep] + 7000, img[keep], err[keep]
    x, y = wcs.world_to_pixel_values(*STAR)
    yy, xx = np.indices(img.shape[1:])
    dist = np.hypot(xx - x, yy - y)
    remote = dist >= 4
    ref = np.median(img, axis=0)
    gy, gx = np.gradient(ref)
    remote &= ref <= np.median(ref[remote])
    design = np.stack(
        [np.ones_like(xx), xx - x, yy - y, gx, gy, ref - ref.mean()], axis=-1
    )
    diff = img - ref
    for j, frame in enumerate(diff):
        use = remote.copy()
        for _ in range(3):
            beta = np.linalg.lstsq(design[use], frame[use], rcond=None)[0]
            residual = frame - design @ beta
            sigma = 1.4826 * np.median(abs(residual[use] - np.median(residual[use])))
            use = remote & (abs(residual) < 3.5 * sigma)
        diff[j] = residual
    nuisance = []
    for ix in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > 0.4) + 1):
        tt = t[ix] - t[ix].mean()
        for k in range(3):
            v = np.zeros(len(t))
            v[ix] = tt**k
            nuisance.append(v)
    results = {}
    for radius in (0.7, 1.2, 1.8, 2.5):
        aperture = dist <= radius
        f = diff[:, aperture].sum(axis=1)
        e = np.sqrt((err[:, aperture] ** 2).sum(axis=1))
        results[str(radius)] = fit(t, f, e, np.stack(nuisance, axis=1))
    return results


def main(cache, output):
    result = {
        "source": "OGLE-LMC-CEP-1429",
        "catalog_period_days": PERIOD,
        "position_deg": STAR,
        "tess": {},
        "ogle": {},
    }
    for sector in (1, 2, 3, 5, 8, 12):
        field = dict(ra_deg=CENTER[0], dec_deg=CENTER[1], size_pixels=13)
        cube, meta = acquire(field, dict(sector=sector), cache)
        result["tess"][str(sector)] = {
            "product": meta,
            "aperture_radius_pixels": tess(cube),
        }
    for campaign, directory in [(3, "ogle3/OIII-CVS"), (4, "ogle4/OCVS")]:
        for band in "IV":
            url = f"https://ogle.astrouw.edu.pl/ogle/{directory}/lmc/cep/phot/{band}/OGLE-LMC-CEP-1429.dat"
            path = cache / f"cep{campaign}_{band}.dat"
            if not path.exists():
                atomic_bytes(path, download(url, 1_000_000))
            raw = path.read_bytes()
            t, mag, err = np.loadtxt(io.BytesIO(raw)).T
            flux = 10 ** (-0.4 * (mag - np.median(mag)))
            result["ogle"][f"{campaign}-{band}"] = {
                "url": url,
                "sha256": hashlib.sha256(raw).hexdigest(),
                **fit(t, flux, flux * np.hypot(err, 0.005)),
            }
    atomic_bytes(
        output, (json.dumps(result, indent=2, allow_nan=False) + "\n").encode()
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    main(a.cache, a.output)
