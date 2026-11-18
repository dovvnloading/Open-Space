"""Check LYACAL pointing identity and empirical scattering subtraction."""

from argparse import ArgumentParser
from hashlib import file_digest, sha256
import io
import json
from pathlib import Path
import tarfile

from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import numpy as np

from check_contrasts import ROWS, bracket, coefficients, line_profiles, measure


TABLE2_MET = 555002780
TABLE2_ICRS = SkyCoord(10.095 * u.deg, -35.029 * u.deg)


def frame(hdus, content, path, calibrated, cols):
    header = hdus[0].header
    record = {
        "path": path, "sha256": sha256(content).hexdigest(), "met": int(header["MET"]),
        "duration": float(header["EXPTIME"]), "door": int(header["APRDRPOS"]),
        "ra": float(header["SPCBRRA"]), "dec": float(header["SPCBRDEC"]),
    }
    record["observation_description"] = header.get("SAPDESC")
    data = hdus["WAVELENGTH_IMAGE"].data if calibrated else hdus[0].data
    if data.shape != (32, 1024):
        raise ValueError(f"Unexpected Alice shape: {path}")
    return record, data[:, cols].copy()


def current_frames(archive_dir, sources, cols):
    pairs = {}
    for kind in ("raw", "calibrated"):
        source = sources[kind]
        path = archive_dir / source["file"]
        with path.open("rb") as stream:
            if file_digest(stream, "sha256").hexdigest() != source["raw_sha256"]:
                raise ValueError(f"Input digest mismatch: {path.name}")
        with tarfile.open(path, "r|gz") as archive:
            for member in archive:
                if not member.isfile() or not member.name.endswith(".fit"):
                    continue
                content = archive.extractfile(member).read()
                with fits.open(io.BytesIO(content), memmap=False) as hdus:
                    header = hdus[0].header
                    if header.get("SAPDESC") != "ALICE CUVB Ly-alpha":
                        continue
                    if header["EXPTIME"] != 3600 and header["MET"] != TABLE2_MET:
                        continue
                    entry = frame(hdus, content, member.name, kind == "calibrated", cols)
                    pairs.setdefault(entry[0]["met"], {})[kind] = entry
    return pairs


def legacy_frames(directory, manifest, cols):
    pairs = {}
    for source in manifest:
        path = directory / source["file"]
        content = path.read_bytes()
        if sha256(content).hexdigest() != source["sha256"]:
            raise ValueError(f"Input digest mismatch: {path.name}")
        sky_calibrated = path.name.endswith("_sci.fit")
        with fits.open(io.BytesIO(content), memmap=False) as hdus:
            record, data = frame(hdus, content, source["file"], sky_calibrated, cols)
        kind = "calibrated" if sky_calibrated else "raw"
        pairs.setdefault(record["met"], {})[kind] = (record, data)
    return pairs


def matching(pairs, door, n):
    selected = []
    for met, pair in sorted(pairs.items()):
        raw, counts = pair["raw"]
        cal, wave = pair["calibrated"]
        if met == TABLE2_MET or raw["door"] != door:
            continue
        if any(raw[key] != cal[key] for key in ("door", "duration", "ra", "dec")):
            raise ValueError(f"Raw/calibrated mismatch: {met}")
        selected.append((dict(raw, calibrated_path=cal["path"], calibrated_sha256=cal["sha256"]), counts, wave))
    if len(selected) != n:
        raise ValueError(f"Unexpected number of frames with door state {door}: {len(selected)}")
    return selected


def direction(coordinate):
    return {"ra": float(coordinate.icrs.spherical.lon.deg), "dec": float(coordinate.icrs.spherical.lat.deg),
            "l": float(coordinate.galactic.l.deg), "b": float(coordinate.galactic.b.deg)}


def mean_boresight(entries):
    coords = SkyCoord([row["ra"] for row, *_ in entries] * u.deg,
                      [row["dec"] for row, *_ in entries] * u.deg)
    xyz = coords.cartesian.xyz.mean(axis=1)
    return SkyCoord(*xyz, representation_type="cartesian", frame="icrs")


def rate(entries):
    counts = np.sum([data for _, data, _ in entries], axis=0, dtype=np.float64)
    duration = sum(row["duration"] for row, *_ in entries)
    return counts / duration, counts / duration**2


def template(old, current, mask, darks=None):
    numerator, variance = rate(old)
    latest, latest_var = rate(current)
    numerator -= latest
    variance += latest_var
    if darks is not None:
        for sign, group in zip((-1, 1), darks):
            dark, dark_var = rate(group)
            numerator += sign * dark
            variance += dark_var
    normalization = float(np.sum(numerator * mask))
    if normalization <= 0:
        raise ValueError("The earlier Ly-alpha signal must exceed the later one")
    return numerator / normalization, normalization, variance


def sensitivity(observations, counts, wave, mask, trials):
    pairs, _, darks = bracket(observations)
    sky = sorted(pairs)
    lya_profiles = np.broadcast_to(mask, wave.shape)
    lya = []
    for index in sky:
        sky_rate = np.sum(mask * counts[index]) / observations[index]["duration"]
        dark_rate = sum(weight * np.sum(mask * counts[i]) / observations[i]["duration"]
                        for i, weight in pairs[index])
        lya.append(sky_rate - dark_rate)
    result = {}
    for line in ("830", "1135", "1135_wide", "1135_curved", "1549"):
        profiles = line_profiles(wave, ROWS["stem"], line)
        mean_sky = np.einsum("i,ijk->jk", lya, profiles[sky]) / len(sky)
        original, _ = coefficients(profiles, observations, pairs, darks, sky)
        row = {"before": measure(counts, original)}
        for name, (pattern, normalization, variance) in trials.items():
            projection = np.einsum("ijk,jk->i", profiles, pattern)
            projected = profiles - projection[:, None, None] * lya_profiles
            local, _ = coefficients(projected, observations, pairs, darks, sky)
            conditioned = measure(counts, local)
            derivative = (mean_sky - np.sum(mean_sky * pattern) * mask) / normalization
            template_error = float(np.sqrt(np.sum(derivative**2 * variance)))
            row[name] = dict(conditioned, template_error=template_error,
                             combined_count_error=float(np.hypot(template_error, conditioned["count_error"])),
                             subtracted_rate=float(np.sum(mean_sky * pattern)))
        result[line] = row
    return result


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("sources", type=Path)
    parser.add_argument("archive_dir", type=Path)
    parser.add_argument("legacy_sources", type=Path)
    parser.add_argument("legacy_dir", type=Path)
    parser.add_argument("science_windows", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    obs = json.loads((args.science_windows / "observations.json").read_bytes())["records"]
    with np.load(args.science_windows / "windows.npz", allow_pickle=False) as study:
        counts, wave, cols = study["counts"], study["wave"], study["columns"]
    modern = current_frames(args.archive_dir, json.loads(args.sources.read_bytes()), cols)
    legacy = legacy_frames(args.legacy_dir, json.loads(args.legacy_sources.read_bytes()), cols)
    old, latest = matching(legacy, 2, 6), matching(modern, 2, 10)
    old_dark, latest_dark = matching(legacy, 1, 3), matching(modern, 1, 11)
    modern_wave = np.median([w for _, _, w in latest], axis=0)
    legacy_wave = np.median([w for _, _, w in old], axis=0)
    mask = (modern_wave >= 1195) & (modern_wave <= 1235)
    mask[list(range(6)) + list(range(16, 32))] = False
    short = modern[TABLE2_MET]["raw"][0]
    short_direction = SkyCoord(short["ra"] * u.deg, short["dec"] * u.deg)
    old_direction, latest_direction = mean_boresight(old), mean_boresight(latest)
    trials = {
        "all_shared_dark": template(old, latest, mask),
        "all_separate_dark": template(old, latest, mask, (old_dark, latest_dark)),
        "early_halves": template(old[:3], latest[:5], mask),
        "late_halves": template(old[3:], latest[5:], mask),
    }
    result = {
        "pointings": {
            "short": dict(short, **{"galactic": direction(short_direction)}),
            "old": direction(old_direction), "new": direction(latest_direction),
            "short_table2_error_arcsec": float(short_direction.separation(TABLE2_ICRS).arcsec),
            "short_new_sep_deg": float(short_direction.separation(latest_direction).deg),
            "old_new_sep_deg": float(old_direction.separation(latest_direction).deg),
        },
        "line_window_old_minus_new_angstrom": float(np.median((legacy_wave-modern_wave)[
            (modern_wave > 1120) & (modern_wave < 1150) & np.isin(np.arange(32)[:, None], ROWS["stem"])])),
        "lya_window_old_minus_new_angstrom": float(np.median((legacy_wave-modern_wave)[mask])),
        "observations": {name: [row for row, *_ in group] for name, group in
                         (("old_sky", old), ("new_sky", latest), ("old_dark", old_dark), ("new_dark", latest_dark))},
        "lya_normalization_rate": {name: value[1] for name, value in trials.items()},
        "contrasts": sensitivity(obs, counts, wave, mask, trials),
    }
    args.output_dir.mkdir()
    with (args.output_dir / "template-results.json").open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
