"""Recover and localize the LMC-4 signal from coordinate-centered TESS images."""

import argparse
import io
import tarfile
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import least_squares

from open_space.ffi_archive import acquire
from open_space.tess_archive import atomic_bytes, download

CENTER = (77.0550166, -69.2370982)
CONTROL = (76.99845833333333, -69.20333333333333)
CONTROL_FREQUENCY = 2 / 2.9590513  # OGLE-LMC-ECL-08126
WIDE_SECTORS = (1, 2, 3, 5, 7, 12, 13, 27, 28)


def offsets(ra, dec):
    return np.array(
        [
            (ra - CENTER[0]) * np.cos(np.deg2rad(CENTER[1])) * 3600,
            (dec - CENTER[1]) * 3600,
        ]
    )


def condition(cube):
    t, img, err, wcs = cube
    # The local aperture and spatial fit cannot use missing detector pixels.
    keep = np.all(np.isfinite(img) & np.isfinite(err) & (err > 0), axis=(1, 2))
    t, img = t[keep], img[keep]
    if len(t) < 300:
        return None
    x, y = wcs.world_to_pixel_values(*CENTER)
    yy, xx = np.indices(img.shape[1:])
    dist = np.hypot(xx - x, yy - y)
    ref = np.median(img, axis=0)
    remote = dist >= 4
    remote &= ref <= np.median(ref[remote])
    gy, gx = np.gradient(ref)
    basis = np.stack(
        [np.ones_like(xx), xx - x, yy - y, gx, gy, ref - ref.mean()], axis=-1
    )
    clean = img - ref
    scatter = []
    for j in range(len(t)):
        selected = remote.copy()
        for _ in range(3):
            beta = np.linalg.lstsq(basis[selected], clean[j][selected], rcond=None)[0]
            residual = clean[j] - basis @ beta
            sigma = 1.4826 * np.median(
                abs(residual[selected] - np.median(residual[selected]))
            )
            selected = remote & (abs(residual) < 3.5 * sigma)
        clean[j] = residual
        scatter.append(np.median(abs(residual[remote])))
    scatter = np.array(scatter)
    typical = np.median(scatter)
    keep = scatter < typical + 6 * 1.4826 * np.median(abs(scatter - typical))
    t, clean = t[keep], clean[keep]
    for ix in np.split(np.arange(len(t)), np.flatnonzero(np.diff(t) > 0.4) + 1):
        if len(ix) < 2:
            continue
        st = t[ix]
        anchors = np.linspace(st[0], st[-1], max(2, int(np.ptp(st) / 0.25) + 1))
        med = np.stack(
            [np.median(clean[ix][abs(st - a) < 1.5], axis=0) for a in anchors]
        )
        left = np.clip(np.searchsorted(anchors, st) - 1, 0, len(anchors) - 2)
        frac = (st - anchors[left]) / (anchors[left + 1] - anchors[left])
        clean[ix] -= (1 - frac)[:, None, None] * med[left] + frac[:, None, None] * med[
            left + 1
        ]
    return t, clean, clean[:, dist < 1.2].sum(axis=1), wcs


def peak(sets):
    time = np.concatenate([v[0] for v in sets])
    amplitude = np.concatenate([v[2] for v in sets])
    frequency, power = LombScargle(time, amplitude).autopower(
        minimum_frequency=0.6, maximum_frequency=1.3, samples_per_peak=12
    )
    i = np.argmax(power)
    return dict(
        frequency_per_day=frequency[i],
        period_days=1 / frequency[i],
        variance_fraction=power[i],
        grid_points=len(frequency),
        baseline_days=np.ptp(time),
    )


def localize(cube, frequency, location, radius=2.55):
    t, images, _, wcs = cube
    phi = 2 * np.pi * frequency * (t - 1500)
    design = np.column_stack([np.ones(len(t)), np.cos(phi), np.sin(phi)])
    coeff = np.linalg.lstsq(design, images.reshape(len(t), -1), rcond=None)[0][
        1:
    ].reshape((2,) + images.shape[1:])
    cx, cy = wcs.world_to_pixel_values(*location)
    yy, xx = np.indices(images.shape[1:])
    if min(cx, cy, xx.max() - cx, yy.max() - cy) < radius:
        return {"status": "edge"}
    dist = np.hypot(xx - cx, yy - cy)
    vector = coeff[:, dist < 1.2].sum(axis=1)
    amplitude = np.linalg.norm(vector)
    projected = np.tensordot(vector / amplitude, coeff, axes=(0, 0))
    use = dist < radius

    def residual(params):
        amp, xc, yc, sx, sy, rho, level = params
        dx = (xx[use] - xc) / sx
        dy = (yy[use] - yc) / sy
        return (
            amp
            * np.exp(-0.5 * (dx * dx + dy * dy - 2 * rho * dx * dy) / (1 - rho * rho))
            + level
            - projected[use]
        )

    scale = max(0.02, np.std(projected[dist > 3]))
    fit = least_squares(
        residual,
        [np.ptp(projected[use]), cx, cy, 0.75, 0.75, 0, np.median(projected[dist > 3])],
        bounds=(
            [0, cx - 2, cy - 2, 0.4, 0.4, -0.6, -np.inf],
            [np.inf, cx + 2, cy + 2, 2.5, 2.5, 0.6, np.inf],
        ),
        loss="soft_l1",
        f_scale=scale,
    )
    sky = wcs.pixel_to_world_values(*fit.x[1:3])
    east, north = offsets(*sky)
    return {
        "status": "measured",
        "ra_deg": float(sky[0]),
        "dec_deg": float(sky[1]),
        "east_arcsec": east,
        "north_arcsec": north,
        "aperture_semiamplitude_e_per_s": amplitude,
        "phase_cycles": np.arctan2(vector[1], vector[0]) / (2 * np.pi),
        "psf_sigma_pixel": fit.x[3:5].tolist(),
        "gaussian_bound_active": bool(np.any(fit.active_mask[:6])),
        "median_abs_pixel_residual_e_per_s": np.median(abs(fit.fun)),
    }


def resolved_photometry(cache, history, frequency):
    data = {}

    def fetch(item, limit):
        info = item["receipt"]
        path = cache / info["sha256"]
        if not path.exists():
            atomic_bytes(path, download(info["url"], limit))
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != info["sha256"]:
            raise ValueError("Resolved series differs from pinned bytes")
        return raw

    def series(t, mag, err):
        t, mag, err = map(np.asarray, (t, mag, err))
        use = np.isfinite(t + mag + err) & (err > 0) & (err < 0.2)
        t, mag, err = t[use], mag[use], err[use]
        measurements = {}
        baseline = np.array([np.median(mag[abs(t - at) < 12]) for at in t])
        for key, values in [("raw", mag), ("running_median_24d", mag - baseline)]:
            ls = LombScargle(t, values, err)
            coeff = ls.model_parameters(frequency)
            measurements[key] = {
                "power": float(ls.power(frequency)),
                "semiamplitude_mag": float(np.hypot(*coeff[1:])),
            }
        return {"epochs": len(t), "span_days": np.ptp(t), "fits": measurements}

    name = "OGLE050813.83-691422.9.mag"
    ogle = history["ground_based_sources"][name]
    raw = fetch(ogle, 50_000_000)
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        member = archive.extractfile(name).read()
    if hashlib.sha256(member).hexdigest() != ogle["member_sha256"]:
        raise ValueError("wrong OGLE-II member")
    data["OGLE-II I"] = dict(
        **series(*np.loadtxt(io.BytesIO(member)).T),
        source=ogle["receipt"],
        member=name,
        member_sha256=ogle["member_sha256"],
    )
    macho = history["ground_based_sources"]["MACHO 1.4532.22"]
    rows = list(csv.DictReader(io.StringIO(fetch(macho, 5_000_000).decode())))
    for band in ("r", "b"):
        selected = [
            r
            for r in rows
            if r[band + "cut"] == "1"
            and r[band + "mag"]
            and r[band + "err"]
            and 0 < float(r[band + "err"]) < 0.08
        ]
        t, m, e = np.array(
            [
                [float(r[k]) for k in ("dateobs", band + "mag", band + "err")]
                for r in selected
            ]
        ).T
        good = abs(m - np.median(m)) < 6 * 1.4826 * np.median(abs(m - np.median(m)))
        data["MACHO " + band] = dict(
            **series(t[good], m[good], e[good]), source=macho["receipt"]
        )
    return data


def analyze(cache, output, reference):
    result = {
        "center_icrs_deg": CENTER,
        "search_frequency_interval_per_day": [0.6, 1.3],
        "train_sectors": list(range(1, 7)),
        "temporal_check_sectors": list(range(7, 14)),
        "control": {
            "object": "OGLE-LMC-ECL-08126",
            "position_deg": CONTROL,
            "frequency_per_day": CONTROL_FREQUENCY,
        },
        "cutouts": [],
    }
    train_frequency = None
    for size, sectors in [(13, range(1, 14)), (21, WIDE_SECTORS)]:
        sets = {}
        for s in sectors:
            spec = {"sector": s}
            if reference:
                prior = next(
                    x
                    for x in reference["cutouts"]
                    if x["size"] == size and x["sector"] == s
                )
                if "product" in prior:
                    spec["science_sha256"] = prior["product"]["science_sha256"]
            field = dict(ra_deg=CENTER[0], dec_deg=CENTER[1], size_pixels=size)
            row = {"size": size, "sector": s}
            result["cutouts"].append(row)
            try:
                cube, meta = acquire(field, spec, cache)
            except ValueError as exc:
                if str(exc) != "Insufficient valid FFI cadences":
                    raise
                row.update(status="insufficient_detector_coverage", reason=str(exc))
                continue
            row["product"] = meta
            measurement = condition(cube)
            if measurement is None:
                row["status"] = "insufficient_complete_image_rows"
                continue
            sets[s] = measurement
            row.update(status="measured", n=len(measurement[0]))
        if size == 13:
            for label, members in [
                ("first_six", range(1, 7)),
                ("remaining_prime", range(7, 14)),
                ("full_prime", range(1, 14)),
            ]:
                result[label] = peak([sets[i] for i in members if i in sets])
            train_frequency = result["first_six"]["frequency_per_day"]
        for rec in result["cutouts"]:
            s = rec["sector"]
            if rec["size"] != size or s not in sets:
                continue
            rec["position"] = localize(sets[s], train_frequency, CENTER)
            if size == 21:
                rec["control_position"] = localize(sets[s], CONTROL_FREQUENCY, CONTROL)
                rec["target_patch_3_2"] = localize(
                    sets[s], train_frequency, CENTER, 3.2
                )
    sources = (
        Path(__file__).resolve().parents[1] / "results" / "ffi-independent-ground.json"
    )
    history = json.loads(sources.read_bytes())
    known = history["simbad"]["objects"]
    result["resolved_at_training_frequency"] = resolved_photometry(
        cache, history, train_frequency
    )
    result["catalog_sources_sha256"] = hashlib.sha256(sources.read_bytes()).hexdigest()
    rows = result["cutouts"]
    primary = [r for r in rows if r["size"] == 13 and r["status"] == "measured"]
    xy = np.array(
        [[r["position"][k] for k in ("east_arcsec", "north_arcsec")] for r in primary]
    )
    center = np.median(xy, axis=0)
    result["thirteen_pixel_summary"] = {
        "median_east_north_arcsec": center.tolist(),
        "east_north_iqr_arcsec": np.percentile(xy, [25, 75], axis=0).tolist(),
        "separation_from_median_arcsec": {
            obj["main_id"]: float(
                np.linalg.norm(center - offsets(float(obj["ra"]), float(obj["dec"])))
            )
            for obj in known
        },
    }
    gaia_url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(
        {
            "-source": "I/355/gaiadr3",
            "-c": "77.0550166 -69.2370982",
            "-c.rs": "30",
            "-out.max": "1000",
            "-out": "Source,RA_ICRS,DE_ICRS,Gmag,VarFlag",
            "-out.add": "_r",
        }
    )
    path = cache / "lmc4_gaia.tsv"
    if not path.exists():
        atomic_bytes(path, download(gaia_url, 1_000_000))
    raw = path.read_bytes()
    lines = [
        l for l in raw.decode().splitlines() if l.strip() and not l.startswith("#")
    ]
    if b"#Name: I/355/gaiadr3" not in raw or not lines[2].startswith("---"):
        raise ValueError("Missing Gaia source rows")
    nearby = []
    for r in csv.DictReader([lines[0], *lines[3:]], delimiter="\t"):
        a, d = map(float, (r["RA_ICRS"], r["DE_ICRS"]))
        nearby.append(
            dict(
                gaia_dr3=r["Source"].strip(),
                ra_deg=a,
                dec_deg=d,
                gmag=float(r["Gmag"]) if r["Gmag"].strip() else None,
                distance_from_median_arcsec=float(
                    np.linalg.norm(center - offsets(a, d))
                ),
            )
        )
    if len(nearby) >= 1000:
        raise ValueError("Truncated Gaia cone")
    result["gaia_cone"] = {
        "url": gaia_url,
        "raw_sha256": hashlib.sha256(raw).hexdigest(),
        "cone_rows": len(nearby),
        "sources_within_5_arcsec": sorted(
            (s for s in nearby if s["distance_from_median_arcsec"] <= 5),
            key=lambda s: s["distance_from_median_arcsec"],
        ),
    }
    atomic_bytes(
        output, (json.dumps(result, indent=2, allow_nan=False) + "\n").encode()
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--verify",
        type=Path,
        help="Evidence JSON whose science digests will be enforced",
    )
    args = parser.parse_args()
    analyze(
        args.cache,
        args.output,
        json.loads(args.verify.read_bytes()) if args.verify else None,
    )
