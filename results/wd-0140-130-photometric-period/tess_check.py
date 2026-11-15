"""Forward-project the ZTF r-band waveform onto independent TESS pixel photometry."""

import argparse
import hashlib
import io
import json
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from analyze import HERE, REF, read, clean_ztf, fit, subset

POSITION = (25.68910759346, 13.26281275445)
SECTORS = (42, 43, 70)


def get_file(cache, sector):
    file = cache / f"tess{sector:04}.fits"
    if not file.exists():
        url = (
            "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=25.68911&dec=13.26281"
            f"&y=7&x=7&units=px&sector={sector}"
        )
        with urllib.request.urlopen(url, timeout=180) as stream:
            with zipfile.ZipFile(io.BytesIO(stream.read())) as archive:
                file.write_bytes(
                    archive.read(
                        next(k for k in archive.namelist() if k.endswith(".fits"))
                    )
                )
    return file


def extract(file):
    with fits.open(file) as hdus:
        w = WCS(hdus[2].header)
        xy = np.array(w.world_to_pixel_values(*POSITION))
        x, y = np.rint(xy).astype(int)
        table = hdus[1].data
        t = table["TIME"].astype(float) + 2457000 - REF
        cube = table["FLUX"].astype(float)
        use = (table["QUALITY"] == 0) & np.isfinite(t) & np.isfinite(cube[:, y, x])
        t, cube = t[use], cube[use]
    cube -= np.nanmedian(cube, axis=0)
    yy, xx = np.indices(cube.shape[1:])
    outer = cube[:, np.hypot(xx - xy[0], yy - xy[1]) > 2.5]
    outer = outer[:, np.all(np.isfinite(outer), axis=0)]
    u = np.linalg.svd(outer, full_matrices=False)[0]
    base = np.column_stack([np.ones(len(t)), t - t.mean(), u[:, :5]])
    central = cube[:, y, x]
    residual = central - base @ np.linalg.lstsq(base, central, rcond=None)[0]
    good = abs(residual - np.median(residual)) < 5 * 1.4826 * np.median(
        abs(residual - np.median(residual))
    )
    neighbors = {
        str(key): list(map(float, w.world_to_pixel_values(ra, dec)))
        for key, ra, dec in [
            ("fainter_G20.37", 25.68217233854, 13.25922135331),
            ("brighter_G15.24", 25.70502081947, 13.26907429914),
        ]
    }
    return t[good], cube[good], base[good], x, y, xy, neighbors


def sector_fit(t, cube, base, x, y, freq, phase):
    flux = cube[:, y, x]
    # Flux is anticorrelated with magnitude, so the predicted phase includes pi.
    template = np.cos(2 * np.pi * freq * t - phase)

    def residual(v):
        return v - base @ np.linalg.lstsq(base, v, rcond=None)[0]

    signal = residual(template)
    nullflux = residual(flux)
    denom = signal @ signal
    numerator = nullflux * signal
    _, inv = np.unique(np.floor(t).astype(int), return_inverse=True)
    by_day = np.bincount(inv, weights=numerator)
    scale = np.std(nullflux, ddof=base.shape[1])
    coords = np.column_stack(
        [base, np.cos(2 * np.pi * freq * t), np.sin(2 * np.pi * freq * t)]
    )
    cc = np.linalg.lstsq(coords, np.nan_to_num(cube.reshape(len(t), -1)), rcond=None)[
        0
    ][-2:]
    image = (np.array([np.cos(phase), np.sin(phase)]) @ cc).reshape(cube.shape[1:])
    peak = np.unravel_index(np.argmax(image), image.shape)
    groups = np.unique(np.floor(t).astype(int))
    block_scores = {}
    for offset in (0, 1):
        _, group = np.unique((groups + offset) // 2, return_inverse=True)
        block_scores[f"two_day_offset{offset}"] = np.bincount(
            group, weights=by_day / scale**2
        )
    block_scores["one_day"] = by_day / scale**2
    return {
        "n": len(t),
        "n_days": len(by_day),
        "fixed_phase_coefficient": float(sum(numerator) / denom),
        "signed_normalized_score": float(sum(numerator) / np.sqrt(denom) / scale),
        "brightest_projected_pixel_xy": [int(peak[1]), int(peak[0])],
    }, block_scores


def run(cache):
    cache.mkdir(parents=True, exist_ok=True)
    info = json.loads((HERE / "measurement.json").read_text())
    d = clean_ztf(read("ztf.csv.gz"))
    r = subset(d, d[3] == "zr")
    tests = {}
    for k in ["signal", "plus_sidereal", "minus_sidereal"]:
        f = info["aliases"][k]["frequency"]
        beta = fit(*r, f)[1]
        tests[k] = (f, np.arctan2(-beta[2], -beta[1]))
    result = {}
    days = {
        k: {method: [] for method in ["one_day", "two_day_offset0", "two_day_offset1"]}
        for k in tests
    }
    for sector in SECTORS:
        file = get_file(cache, sector)
        t, cube, base, x, y, xy, neighbors = extract(file)
        fits_hash = hashlib.sha256(file.read_bytes()).hexdigest()
        result[str(sector)] = {
            "target_xy": xy.tolist(),
            "central_pixel_xy": [int(x), int(y)],
            "start_btjd": float(min(t) + 1000),
            "end_btjd": float(max(t) + 1000),
            "median_cadence_seconds": round(float(np.median(np.diff(t)) * 86400)),
            "cutout_sha256": fits_hash,
            "neighbor_locations_xy": neighbors,
        }
        for k, (f, phase) in tests.items():
            data, blocks = sector_fit(t, cube, base, x, y, f, phase)
            result[str(sector)][k] = data
            if k == "signal":
                # The approximate inter-survey time-scale difference is smaller than
                # 90 s. Recalculate both endpoints without tuning the waveform.
                data["clock_sensitivity_90s"] = [
                    sector_fit(
                        t, cube, base, x, y, f, phase + sgn * 2 * np.pi * f * 90 / 86400
                    )[0]["fixed_phase_coefficient"]
                    for sgn in [-1, 1]
                ]
                control = {}
                for neighbor, coords in neighbors.items():
                    xx, yy = np.rint(coords).astype(int)
                    if (
                        0 <= xx < cube.shape[2]
                        and 0 <= yy < cube.shape[1]
                        and np.all(np.isfinite(cube[:, yy, xx]))
                    ):
                        control[neighbor] = sector_fit(t, cube, base, xx, yy, f, phase)[
                            0
                        ]["fixed_phase_coefficient"]
                data["neighbor_projected_coefficients"] = control
            for method, value in blocks.items():
                days[k][method].extend(value.tolist())
    rng = np.random.default_rng(120407)
    stats = {}
    for k, methods in days.items():
        stats[k] = {}
        for method, blocks in methods.items():
            scores = np.array(blocks)
            hits = 0
            for _ in range(200):
                hits += int(
                    np.sum(
                        (rng.choice([-1, 1], (500, len(scores))) @ scores)
                        >= sum(scores)
                    )
                )
            stats[k][method] = {
                "blocks": len(scores),
                "trials": 100000,
                "as_or_more_extreme": hits,
                "one_sided_p": (hits + 1) / 100001,
            }
    return {"tess_sectors": result, "night_block_fixed_phase_tests": stats}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cache", type=Path)
    print(json.dumps(run(parser.parse_args().cache), indent=2))
