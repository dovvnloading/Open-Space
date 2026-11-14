"""Replay the four-visit emission-line screen and its targeted controls."""

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import angular_separation
from astropy.io import fits
from astropy.table import Table
from scipy.ndimage import convolve1d, median_filter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
VISITS = ((6497, 56329), (6498, 56565), (6503, 56564), (6509, 56486))
LINES = np.array([5008.239, 4960.295, 4862.683, 3728.48])
C = 299792.458
SKY_COORD = (338.5431068618596, 32.4453897077225)


def file_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_inputs(cache):
    cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / "data/sdss-superposition-inputs.json").read_text())
    for name, receipt in manifest.items():
        dest = cache / name
        if not dest.exists():
            if name == "desi.csv":
                shutil.copyfile(ROOT / "data/sdss-superposition-desi.csv", dest)
            else:
                with NamedTemporaryFile(dir=cache, delete=False) as tmp:
                    staged = Path(tmp.name)
                try:
                    with (
                        urlopen(receipt["url"], timeout=120) as response,
                        staged.open("wb") as output,
                    ):
                        shutil.copyfileobj(response, output)
                    if file_hash(staged) != receipt["sha256"]:
                        raise ValueError(f"Download hash mismatch for {name}")
                    staged.replace(dest)
                finally:
                    staged.unlink(missing_ok=True)
        if (
            dest.stat().st_size != receipt["size_bytes"]
            or file_hash(dest) != receipt["sha256"]
        ):
            raise ValueError(f"Cache hash mismatch for {name}")


def separation(ra, dec):
    return (
        np.degrees(
            angular_separation(*np.radians(SKY_COORD), np.radians(ra), np.radians(dec))
        )
        * 3600
    )


def search(flux, model, weight, and_mask, wave, table, selected):
    kernel = np.exp(-0.5 * (np.arange(-6, 7) / 1.4) ** 2)
    positions = np.arange(len(wave))
    hits = []
    for i in np.flatnonzero(selected):
        residual = flux[i].astype(float) - model[i]
        ivar = weight[i].astype(float)
        ivar[(and_mask[i] != 0) | ~np.isfinite(residual) | ~np.isfinite(ivar)] = 0
        residual = np.nan_to_num(residual)
        residual -= median_filter(residual, size=41)
        norm = convolve1d(ivar, kernel**2, mode="constant")
        height = convolve1d(ivar * residual, kernel, mode="constant") / np.maximum(
            norm, 1e-12
        )
        score = height * np.sqrt(norm)
        for peak in find_peaks(score, height=4.8, distance=5)[0]:
            redshift = wave[peak] / LINES[0] - 1
            if not table["Z"][i] + 0.03 <= redshift <= 1:
                continue
            pixels = peak + np.log10(LINES / LINES[0]) / 0.0001
            if pixels.min() < 7 or np.any(ivar[np.rint(pixels).astype(int)] <= 0):
                continue
            scores = np.interp(pixels, positions, score)
            if scores[1] >= 2.7 and max(scores[2:]) >= 3.2:
                hits.append(
                    {
                        "fiber": int(i + 1),
                        "z_foreground": float(table["Z"][i]),
                        "z_second": redshift,
                        "scores": scores.tolist(),
                        "fitted_peak_ratio": float(
                            height[peak] / max(height[round(pixels[1])], 0.01)
                        ),
                        "ra": table["PLUG_RA"][i],
                        "dec": table["PLUG_DEC"][i],
                    }
                )
    return hits


def measure_window(wave, flux, ivar, good, center):
    """Weighted Gaussian + local linear continuum. Integrated model flux and diagonal error."""
    take = good & (ivar > 0) & (abs(wave - center) < 15)
    x = wave[take] - center
    if len(x) < 7 or np.sum(abs(x) < 4) < 2:
        return None
    sigma = center * 83 / C
    profile = np.exp(-0.5 * (x / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)
    design = np.column_stack([profile, np.ones_like(x), x / 15])
    weight = ivar[take]
    covariance = np.linalg.inv((design.T * weight) @ design)
    coeff = covariance @ (design.T @ (weight * flux[take]))
    err = np.sqrt(covariance[0, 0])
    return {"flux": coeff[0], "error": err, "snr": coeff[0] / err, "pixels": len(x)}


def oxygen_profile(wave, flux, ivar, good, center):
    """Post-selection observed line width without deconvolving the DESI response."""
    use = good & (ivar > 0) & (abs(wave - center) < 14)
    x, y = wave[use] - center, flux[use]

    def profile(x, integral, shift, width, intercept, slope):
        return (
            integral
            * np.exp(-0.5 * ((x - shift) / width) ** 2)
            / (np.sqrt(2 * np.pi) * width)
            + intercept
            + slope * x
        )

    values, covariance = curve_fit(
        profile,
        x,
        y,
        p0=[30, 0, 1.6, np.median(y), 0],
        sigma=1 / np.sqrt(ivar[use]),
        absolute_sigma=True,
        bounds=([0, -5, 0.5, -np.inf, -np.inf], [np.inf, 5, 6, np.inf, np.inf]),
    )
    error = np.sqrt(np.diag(covariance))
    return {
        "integrated_flux": values[0],
        "integrated_flux_error": error[0],
        "center_offset_A": values[1],
        "center_offset_error_A": error[1],
        "observed_sigma_A": values[2],
        "observed_sigma_error_A": error[2],
    }


def common_mode(wave, flux, model, ivar, mask, table, objtype, plate):
    result = {}
    for center in [8999.07, 9028.07, 9228.74]:
        core = abs(wave - center) < 2.3
        blue = (wave > center - 21) & (wave < center - 10)
        red = (wave > center + 10) & (wave < center + 21)
        subset = core | blue | red
        valid = (ivar[:, subset] > 0) & (mask[:, subset] == 0)
        covered = (
            valid[:, core[subset]].all(axis=1)
            & (valid[:, blue[subset]].sum(axis=1) >= 2)
            & (valid[:, red[subset]].sum(axis=1) >= 2)
        )
        fid = table["FIBERID"]
        galaxy = (
            (table["CLASS"] == "GALAXY") & (table["ZWARNING"] == 0) & (objtype != "SKY")
        )
        # Withhold the three reported fibers only on the discovery visit.
        control = galaxy & (~np.isin(fid, [128, 311, 813]) if plate == 6509 else True)
        index = {}
        for label, rows in [
            ("galaxies", control),
            ("blank_sky", objtype == "SKY"),
            ("standards", objtype == "SPECTROPHOTO_STD"),
            ("reported", galaxy & np.isin(fid, [128, 311, 813]) & (plate == 6509)),
        ]:
            idx = np.flatnonzero(rows & covered)
            if not len(idx):
                continue

            def window(data, band, subset=subset, valid=valid, idx=idx):
                use = band[subset]
                return np.nanmedian(
                    np.where(valid[idx][:, use], data[idx][:, band], np.nan), axis=1
                )

            residual = flux if label == "blank_sky" else flux - model
            amplitude = (
                window(residual, core)
                - (window(residual, blue) + window(residual, red)) / 2
            )
            continuum = (window(model, blue) + window(model, red)) / 2
            bright = continuum > 5
            stats = {"n": len(idx), "median_amplitude": float(np.median(amplitude))}
            if label == "galaxies":
                stats.update(
                    rho=float(spearmanr(amplitude, continuum).statistic),
                    n_bright=int(bright.sum()),
                    bright_positive=int(sum(amplitude[bright] > 0)),
                    bright_fractional_median=float(
                        np.median(amplitude[bright] / continuum[bright])
                    ),
                )
            if label == "reported":
                stats["by_fiber"] = {
                    str(fid[j]): float(val)
                    for j, val in zip(idx, amplitude / continuum)
                }
            index[label] = stats
        result[str(center)] = index
    return result


def plate_results(cache):
    visits = {}
    hits = []
    for plate, mjd in VISITS:
        ppath = cache / ("plate.fits" if plate == 6509 else f"p{plate}.fits")
        with fits.open(ppath) as p, fits.open(cache / f"z{plate}.fits") as z:
            header = p[0].header
            flux, ivar, mask = (p[h].data for h in [0, 1, 2])
            wave = 10 ** (
                header["COEFF0"] + header["COEFF1"] * np.arange(flux.shape[1])
            )
            model, table = z[2].data, z[1].data
            assert header["PLATEID"] == plate and header["MJD"] == mjd
            assert np.array_equal(table["FIBERID"], np.arange(1, 1 + len(table)))
            objtype = p["PLUGMAP"].data["OBJTYPE"]
            selected = (
                (table["CLASS"] == "GALAXY")
                & (table["ZWARNING"] == 0)
                & (table["Z"] > 0.01)
                & (objtype != "SKY")
            )
            for hit in search(flux, model, ivar, mask, wave, table, selected):
                hits.append(dict(plate=plate, mjd=mjd, **hit))
            visits[str(plate)] = {
                "mjd": mjd,
                "selected": int(selected.sum()),
                "common_mode": common_mode(
                    wave, flux, model, ivar, mask, table, objtype, plate
                ),
            }
            if plate == 6509:
                centers = LINES * (
                    1
                    + next(
                        h["z_second"]
                        for h in hits
                        if h["plate"] == 6509 and h["fiber"] == 788
                    )
                )
                neighbor = {}
                for fiber in (787, 788, 789):
                    idx = fiber - 1
                    neighbor[str(fiber)] = [
                        measure_window(
                            wave,
                            flux[idx] - model[idx],
                            ivar[idx],
                            mask[idx] == 0,
                            center,
                        )
                        for center in centers
                    ]
                visits[str(plate)]["adjacent_fibers"] = neighbor
    return visits, hits


def aperture_results(cache, hit):
    centers = LINES * (1 + hit["z_second"])
    outputs = {}
    plot = {}
    with fits.open(cache / "s788.fits") as spec:
        source = spec["SPALL"].data[0]
        outputs["metadata"] = {
            "velocity_dispersion": float(source["VDISP"]),
            "velocity_dispersion_error": float(source["VDISP_ERR"]),
            "seeing50_arcsec": [
                e.header["SEEING50"] for e in spec[4:] if e.name.startswith("R")
            ],
        }
        for h in [spec["COADD"], *[e for e in spec[4:] if e.name.startswith("R")]]:
            d = h.data
            wave = 10 ** d["loglam"].astype(float)
            flux = d["flux"].astype(float)
            good = d["and_mask" if h.name == "COADD" else "mask"] == 0
            if h.name == "COADD":
                flux -= d["model"]
                plot["SDSS foreground aperture"] = (wave, flux, d["ivar"], good)
            outputs[h.name] = [
                measure_window(wave, flux, d["ivar"], good, c) for c in centers
            ]
    redshifts = list(csv.DictReader((cache / "desi.csv").open()))
    for entry in redshifts:
        target = int(entry["targetid"])
        with fits.open(cache / f"desi-{entry['program']}.fits") as f:
            rows = np.flatnonzero(f["FIBERMAP"].data["TARGETID"] == target)
            assert len(rows) == 1
            row = rows[0]
            entry["arcsec_from_sdss"] = float(
                separation(
                    float(entry["mean_fiber_ra"]), float(entry["mean_fiber_dec"])
                )
            )
            exp = f["EXP_FIBERMAP"].data
            entry["exposures"] = [
                {k.lower(): int(row[k]) for k in ("NIGHT", "EXPID", "TILEID", "FIBER")}
                for row in exp[exp["TARGETID"] == target]
            ]
            measurements, segments = [], []
            entry["free_oxygen_profiles"] = []
            for line in centers:
                arm = "R" if line < 7500 else "Z"
                wave = f[f"{arm}_WAVELENGTH"].data
                flux = f[f"{arm}_FLUX"].data[row]
                ivar = f[f"{arm}_IVAR"].data[row]
                good = f[f"{arm}_MASK"].data[row] == 0
                measurements.append(measure_window(wave, flux, ivar, good, line))
                if line in centers[:2]:
                    entry["free_oxygen_profiles"].append(
                        oxygen_profile(wave, flux, ivar, good, line)
                    )
                cut = abs(wave - line) < 30
                segments.append((wave[cut], flux[cut], ivar[cut], good[cut]))
            entry["at_sdss_second_redshift"] = measurements
            label = (
                "DESI offset aperture"
                if entry["program"] == "dark"
                else "DESI foreground aperture"
            )
            plot[label] = tuple(
                np.concatenate([segment[j] for segment in segments]) for j in range(4)
            )
    for fiber in [128, 813]:
        with fits.open(cache / f"s{fiber}.fits") as spec:
            coadd = spec["COADD"].data
            model_wave = 10 ** coadd["loglam"].astype(float)
            positive_excess = {}
            for h in spec[4:]:
                if not h.name.startswith("R"):
                    continue
                d = h.data
                wave = 10 ** d["loglam"].astype(float)
                residual = d["flux"] - np.interp(wave, model_wave, coadd["model"])
                good = (d["ivar"] > 0) & (d["mask"] == 0)
                distance = abs(wave - 9228.74)
                core = good & (distance < 3)
                sides = good & (distance > 10) & (distance < 22)
                if core.sum() >= 2 and sides.sum() >= 4:
                    positive_excess[h.name] = float(
                        np.mean(residual[core]) - np.median(residual[sides])
                    )
            outputs[f"earlier_residual_fiber_{fiber}"] = positive_excess
    return outputs, redshifts, plot


def catalog_matches(cache):
    result = {}
    for name, ra, dec in [
        ("silo.fits", "RA", "DEC"),
        ("desilenses.fits", "TARGET_RA", "TARGET_DEC"),
        ("desi-anomaly.fits", "TARGET_RA", "TARGET_DEC"),
        ("targeted.fits", "TARGET_RA", "TARGET_DEC"),
    ]:
        with fits.open(cache / name) as hdus:
            data = hdus[1].data
            sep = separation(data[ra].astype(float), data[dec].astype(float))
            result[name] = {"rows": len(data), "nearest_arcsec": float(np.nanmin(sep))}
    site_pairs = []
    for gid in ("1462846720", "1286257049", "754966998", "1449192048"):
        with (cache / f"hsu-{gid}.csv").open() as stream:
            site_pairs.extend(
                (float(r["RAdeg"]), float(r["DEdeg"]))
                for r in csv.DictReader(stream)
                if r["RAdeg"].strip()
            )
    result["Hsu author tables"] = {
        "rows": len(site_pairs),
        "nearest_arcsec": float(np.min(separation(*np.array(site_pairs).T))),
    }
    pairs = Table.read(cache / "fullcat.mrt", format="ascii.cds")
    result["fullcat.mrt"] = {
        "rows": len(pairs),
        "nearest_arcsec": float(np.min(separation(pairs["RAdeg"], pairs["DEdeg"]))),
    }
    return result


def emit_plot_data(path, plot, centers):
    with path.open("w", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(
            [
                "aperture",
                "line_rest_A",
                "wavelength_A",
                "flux_density",
                "ivar",
                "unflagged",
            ]
        )
        for label, (wave, flux, ivar, good) in plot.items():
            for rest, center in zip(LINES, centers):
                for j in np.flatnonzero(abs(wave - center) < 26):
                    writer.writerow(
                        [
                            label,
                            rest,
                            f"{wave[j]:.6f}",
                            f"{flux[j]:.7g}",
                            f"{ivar[j]:.7g}",
                            int(good[j]),
                        ]
                    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare_inputs(args.cache)
    args.output.mkdir(parents=True, exist_ok=True)
    visits, hits = plate_results(args.cache)
    target = next(h for h in hits if (h["plate"], h["fiber"]) == (6509, 788))
    sdss, desi, plot = aperture_results(args.cache, target)
    result = {
        "line_order": ["[O III] 5008", "[O III] 4960", "H beta", "[O II] blend"],
        "visits": visits,
        "hits": hits,
        "sdss_forced": sdss,
        "desi_followup": desi,
        "catalog_checks": catalog_matches(args.cache),
    }
    (args.output / "sdss-superposition.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    emit_plot_data(
        args.output / "sdss-superposition-windows.csv",
        plot,
        LINES * (1 + target["z_second"]),
    )


if __name__ == "__main__":
    main()
