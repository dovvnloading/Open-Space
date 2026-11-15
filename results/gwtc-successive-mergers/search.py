"""Reproduce the GWTC-5.0 within-release successive-merger search."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy_healpix import HEALPix

MODELS = ("IMRPhenomXPHM_SpinTaylor", "SEOBNRv5PHM")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def sky_grid(table, sample_pixels, oversample):
    # Decode LVK's multi-order NESTED pixel indices, then query the density
    # at equal-solid-angle subpixel centers. All coordinates are equatorial.
    uniq = table["UNIQ"].astype(np.int64)
    order = np.floor(np.log2(uniq) / 2).astype(int) - 1
    pix = uniq - 4 * 4**order
    high_order = max(order)
    starts = pix << (2 * (high_order - order))
    idx = np.argsort(starts)
    highpix = sample_pixels >> (2 * (15 - high_order))
    cell = idx[np.searchsorted(starts[idx], highpix, side="right") - 1]
    density = table["PROBDENSITY"][cell]
    nlat, nlon = np.array(density.shape) // oversample
    sky = density.reshape(nlat, oversample, nlon, oversample).mean(axis=(1, 3))
    integral = float(sky.mean() * 4 * np.pi)
    return sky / sky.mean(), integral


def load_maps(archive, names, model, nlon, nlat, oversample):
    lon = 2 * np.pi * (np.arange(nlon * oversample) + 0.5) / (nlon * oversample)
    sin_dec = -1 + 2 * (np.arange(nlat * oversample) + 0.5) / (nlat * oversample)
    lon, dec = np.meshgrid(lon, np.arcsin(sin_dec))
    grid = HEALPix(nside=2**15, order="nested")
    sample_pixels = grid.lonlat_to_healpix(lon * u.rad, dec * u.rad)
    maps, integrals = {}, []
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            match = re.search(r"-(GW\d{6}(?:_\d{6})?)-" + model + r"_Skymap_", member.name)
            if not match or match[1] not in names:
                continue
            data = tar.extractfile(member).read()
            if data.startswith(b"\x1f\x8b"):
                data = gzip.decompress(data)
            with fits.open(io.BytesIO(data)) as f:
                assert f[1].header["ORDERING"] == "NUNIQ"
                assert f[1].header["COORDSYS"] == "C"
                sky, integral = sky_grid(f[1].data, sample_pixels, oversample)
            maps[match[1]] = sky
            integrals.append(integral)
    return maps, [min(integrals), max(integrals)]


def discrepancy(event_a, field_a, event_b, field_b):
    a, b = event_a[field_a], event_b[field_b]
    # The archive reports signed offsets enclosing 90% of the posterior.
    sa, sb = ("_lower", "_upper") if a > b else ("_upper", "_lower")
    return 1.644854 * abs(a - b) / np.hypot(event_a[field_a + sa], event_b[field_b + sb])


def rotated_overlaps(skies):
    # Cyclic RA cross-correlation; units give an overlap of one for two
    # uniform sky densities. Declination and the map shapes stay fixed.
    freq = np.fft.rfft(skies, axis=2).transpose(2, 0, 1)
    product = freq @ freq.conj().transpose(0, 2, 1)
    return np.fft.irfft(product, n=skies.shape[2], axis=0) / (skies.shape[1] * skies.shape[2])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("catalog", type=Path)
    p.add_argument("archive", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--model", choices=MODELS, default=MODELS[0])
    p.add_argument("--longitude-bins", type=int, default=256)
    p.add_argument("--latitude-bins", type=int, default=128)
    p.add_argument("--oversample", type=int, default=3)
    p.add_argument("--trials", type=int, default=20000)
    args = p.parse_args()
    catalog = json.loads(args.catalog.read_text())["events"].values()
    fields = ("mass_1_source", "mass_2_source", "final_mass_source", "luminosity_distance")
    eligible = {e["commonName"]: e for e in catalog if all(
        e.get(k + s) is not None for k in fields for s in ("", "_lower", "_upper"))}
    maps, integrals = load_maps(args.archive, eligible, args.model, args.longitude_bins,
                                args.latitude_bins, args.oversample)
    names = sorted(maps, key=lambda n: eligible[n]["GPS"])
    skies = np.array([maps[n] for n in names])
    corr = rotated_overlaps(skies)
    a, b = np.triu_indices(len(names), k=1)
    evidence, weights = [], []
    for i, j in zip(a, b):
        x, y = eligible[names[i]], eligible[names[j]]
        mass_z = [discrepancy(x, "final_mass_source", y, k) for k in fields[:2]]
        component = int(np.argmin(mass_z))
        dist_z = discrepancy(x, "luminosity_distance", y, "luminosity_distance")
        weight = np.exp(-0.5 * (min(mass_z)**2 + dist_z**2))
        weights.append(weight)
        evidence.append(dict(earlier=names[i], later=names[j], component=component + 1,
                             days=(y["GPS"] - x["GPS"]) / 86400,
                             remnant_msun=x[fields[2]], component_msun=y[fields[component]],
                             earlier_mpc=x[fields[3]], later_mpc=y[fields[3]],
                             mass_discrepancy=min(mass_z), distance_discrepancy=dist_z,
                             sky_overlap=corr[0, i, j], score=corr[0, i, j] * weight))
    weights = np.asarray(weights)
    observed = np.max(corr[0, a, b] * weights)
    rng = np.random.default_rng(20261115)
    maxima = []
    for start in range(0, args.trials, 80):
        count = min(80, args.trials - start)
        offsets = rng.integers(0, args.longitude_bins, size=(count, len(names)))
        shifts = (offsets[:, b] - offsets[:, a]) % args.longitude_bins
        maxima.extend((corr[shifts, a, b] * weights).max(axis=1))
    n_exceed = int(np.count_nonzero(maxima >= observed))
    summary = dict(model=args.model, catalog_rows=len(catalog), has_parameters=len(eligible),
                   with_sky_maps=len(names), pairs=len(a), raster_integral_range=integrals,
                   longitude_bins=args.longitude_bins, latitude_bins=args.latitude_bins,
                   oversample=args.oversample, seed=20261115, trials=args.trials,
                   max_score=observed, null_exceedances=n_exceed,
                   global_p=(n_exceed + 1) / (args.trials + 1),
                   null_quantiles_50_90_99=np.quantile(maxima, [0.5, 0.9, 0.99]).tolist(),
                   sha256={p.name: digest(p) for p in (args.catalog, args.archive)})
    args.output.mkdir(exist_ok=True, parents=True)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    evidence.sort(key=lambda e: -e["score"])
    with (args.output / "top_pairs.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(evidence[0]), lineterminator="\n")
        writer.writeheader()
        for row in evidence[:20]:
            writer.writerow({k: round(v, 6) if isinstance(v, float) else v for k, v in row.items()})
    print(json.dumps({k: summary[k] for k in ("pairs", "global_p", "max_score")}))
    print(evidence[0]["earlier"], evidence[0]["later"])


if __name__ == "__main__":
    main()
