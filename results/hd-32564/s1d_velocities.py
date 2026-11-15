"""Measure relative ESPRESSO velocities directly on the barycentric S1D flux grid."""
import argparse
import csv
import hashlib
import io
import json
import tempfile
import time
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from astropy.io.fits import Header

from replay import EPOCH, PREDICTION_HASH, ROOT, design, remove

C = 299792458.0
# Vacuum Angstrom. Exclude the red optical region with stronger telluric bands.
BANDS = [(4200, 4600), (4600, 5000), (5000, 5500)]
PREFIX_BYTES = 13_000_001
BASE_URL = "https://dataportal.eso.org/dataPortal/file/"


def prefix(product, cache):
    path = cache / (product + ".prefix")
    if path.exists():
        data = path.read_bytes()
    else:
        request = Request(BASE_URL + product, headers={"Range": f"bytes=0-{PREFIX_BYTES - 1}"})
        for attempt in range(3):
            try:
                with urlopen(request, timeout=90) as response:
                    data = response.read(PREFIX_BYTES)
                if len(data) != PREFIX_BYTES:
                    raise IOError(f"Short response for {product}")
                break
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(1 + attempt)
        path.write_bytes(data)
    if len(data) != PREFIX_BYTES:
        raise ValueError(f"Incomplete cache for {product}")
    return data


def fits_header(data, start):
    end = next(i for i in range(start, len(data) - 79, 80)
               if data[i:i + 8] == b"END     ") + 80
    raw = data[start:end]
    return Header.fromstring(raw.decode("ascii"), sep=""), (end + 2879) // 2880 * 2880, raw


def spectrum(row, cache):
    data = prefix(row["product"], cache)
    primary, start, raw = fits_header(data, 0)
    table, start, _ = fits_header(data, start)
    if hashlib.sha256(raw).hexdigest() != row["fits_header_sha256"]:
        raise ValueError(f'Primary header changed: {row["product"]}')
    if (primary["OBJECT"] != "HD32564" or primary["SPECSYS"] != "BARYCENT"
            or primary["ESO PRO CATG"] != "S1D_FINAL_A" or table["NAXIS2"] != 1
            or table["TUNIT1"].lower() != "angstrom"):
        raise ValueError("Unexpected product or spectral frame")
    n = int(table["TFORM1"][:-1])
    columns = []
    for i, (name, code, dtype) in enumerate(
            [("WAVE", "D", ">f8"), ("FLUX", "D", ">f8"), ("ERR", "D", ">f8"), ("QUAL", "J", ">i4")], 1):
        if table[f"TTYPE{i}"] != name or table[f"TFORM{i}"] != f"{n}{code}":
            raise ValueError("Unexpected binary table layout")
        column = np.frombuffer(data, dtype=dtype, count=n, offset=start)
        columns.append(column)
        start += column.nbytes
    return (*columns, hashlib.sha256(data).hexdigest())


def order_shifts(flux, err, dv):
    """Fit individual chunk shifts, projecting out a scale and a color term."""
    norm = np.median(flux, axis=1, keepdims=True)
    f, e = flux / norm, err / norm
    template = np.median(f, axis=0)
    grad = -np.gradient(template) / dv
    coordinate = np.linspace(-1, 1, len(template))
    model = np.column_stack([template, coordinate * template, grad])
    weights, shifts = [], []
    for y, sigma in zip(f, e):
        good = np.ones(len(template), dtype=bool)
        for step in range(2):
            weighted = model[good] / sigma[good, None]
            coef = np.linalg.lstsq(weighted, y[good] / sigma[good], rcond=None)[0]
            if step == 0:
                residual = (y - model @ coef) / sigma
                good = np.abs(residual) < 6
                if np.count_nonzero(good) < 0.8 * len(template):
                    return None
        shifts.append(coef[2] / coef[0])
        covariance = np.linalg.inv(weighted.T @ weighted)
        weights.append(coef[0]**2 / covariance[2, 2])
    return np.asarray(shifts), np.asarray(weights)


def measurements(rows, cache, injected=None):
    per_exposure = [spectrum(row, cache) for row in rows]
    wave = per_exposure[0][0]
    for item in per_exposure[1:]:
        np.testing.assert_allclose(item[0], wave, rtol=0, atol=1e-9)
    dv = np.median(np.diff(np.log(wave))) * C
    np.testing.assert_allclose(np.diff(np.log(wave)) * C, dv, rtol=1e-6)
    bands = []
    for lo, hi in BANDS:
        pix = np.flatnonzero((wave >= lo) & (wave < hi))
        f = np.stack([s[1][pix] for s in per_exposure])
        e = np.stack([s[2][pix] for s in per_exposure])
        valid = np.all(np.stack([s[3][pix] == 0 for s in per_exposure]) &
                       np.isfinite(f) & (f > 0) & np.isfinite(e) & (e > 0), axis=0)
        shifts, weights = [], []
        total_pixels = 0
        for start in range(0, len(pix), 512):
            stop = min(start + 512, len(pix))
            if stop - start < 128 or not np.all(valid[start:stop]):
                continue
            chunk = f[:, start:stop]
            if injected is not None:
                x = np.arange(stop - start)
                chunk = np.stack([np.interp(x - velocity / dv, x, y)
                                  for y, velocity in zip(chunk, injected)])
            fit = order_shifts(chunk, e[:, start:stop], dv)
            if fit is None:
                continue
            shift, weight = fit
            shifts.append(shift - np.mean(shift))
            # Identical chunk weights for every visit prevent time-variable weighting.
            weights.append(np.median(weight))
            total_pixels += stop - start
        shifts = np.stack(shifts)
        bands.append((lo, hi, np.average(shifts, axis=0, weights=weights), len(shifts), total_pixels))
    hashes = [s[-1] for s in per_exposure]
    return bands, hashes


def format_csv(rows, bands, hashes):
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(["product", "bjd", "prefix_sha256"] +
                    [f"rv_{lo}_{hi}_m_s" for lo, hi, *_ in bands])
    for i, row in enumerate(rows):
        writer.writerow([row["product"], row["bjd"], hashes[i]] +
                        [f"{band[2][i]:.6f}" for band in bands])
    return buffer.getvalue()


def diagnostics(rows, channels):
    frozen = (ROOT / "harps_predictions.json").read_bytes()
    if hashlib.sha256(frozen).hexdigest() != PREDICTION_HASH:
        raise ValueError("Frozen HARPS file changed")
    models = json.loads(frozen)
    best = min(models, key=lambda m: m["rms"])
    t = np.array([float(r["bjd"]) for r in rows]) - EPOCH
    base = np.column_stack([np.ones(len(t)), (t - t.mean()) / 1000])
    ccf = np.array([float(r["rv"]) * 1000 for r in rows])
    terms = (design(t, best["P"]) * best["params"][:8]).reshape(-1, 4, 2).sum(axis=2)
    rms = lambda a: float(np.sqrt(np.mean(a**2)))
    result = {}
    for name, y in channels.items():
        residual = remove(y, base)
        scores = []
        for m in models:
            predicted = design(t, m["P"]) @ m["params"][:8]
            scores.append((m["P"][0], rms(remove(y - predicted, base))))
        result[name] = {
            "n": len(y), "null_rms": rms(residual),
            "short_only_rms": rms(remove(y - terms[:, 0], base)),
            "full_rms": rms(remove(y - terms.sum(axis=1), base)),
            "without_second_rms": rms(remove(y - terms[:, [0, 2, 3]].sum(axis=1), base)),
            "best_other_alias_rms": min(score for p, score in scores if p != best["P"][0]),
            "ccf_r": float(np.corrcoef(residual, remove(ccf, base))[0, 1]),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", action="store_true", help="download and remeasure all spectral pixels")
    parser.add_argument("--inject", action="store_true", help="repeat extraction with alternating +/-5 m/s shifts")
    parser.add_argument("--cache", type=Path, help="directory to reuse 325 MB of FITS prefixes")
    parser.add_argument("--output", type=Path, help="write remeasured CSV here for comparison")
    args = parser.parse_args()
    if (args.output or args.inject) and not args.extract:
        parser.error("--output and --inject require --extract")
    rows = list(csv.DictReader((ROOT / "espresso.csv").read_text().splitlines()))
    released = ROOT / "s1d-velocities.csv"
    if args.extract:
        with tempfile.TemporaryDirectory() as tmp:
            cache = args.cache or Path(tmp)
            cache.mkdir(exist_ok=True, parents=True)
            bands, hashes = measurements(rows, cache)
            csv_text = format_csv(rows, bands, hashes)
            if args.inject:
                injected = np.where(np.arange(len(rows)) % 2 == 0, 5., -5.)
                shifted, _ = measurements(rows, cache, injected)
                centered = injected - injected.mean()
                for original, perturbed in zip(bands, shifted):
                    gain = centered @ (perturbed[2] - original[2]) / (centered @ centered)
                    print("injected velocity recovery", original[:2], round(float(gain), 5))
        print("band / chunks / eligible pixels", [(b[:2], b[3], b[4]) for b in bands])
        if args.output:
            args.output.write_text(csv_text)
        elif csv_text != released.read_text():
            raise ValueError("Remeasured CSV differs from release")
    else:
        csv_text = released.read_text()
    rv_rows = list(csv.DictReader(csv_text.splitlines()))
    if [(r["product"], r["bjd"]) for r in rows] != [(r["product"], r["bjd"]) for r in rv_rows]:
        raise ValueError("Velocity table does not match exposure manifest")
    channels = {f"{lo}-{hi}": np.array([float(r[f"rv_{lo}_{hi}_m_s"]) for r in rv_rows])
                for lo, hi in BANDS}
    print(json.dumps(diagnostics(rows, channels), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
