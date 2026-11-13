"""Recreate the metadata-only holdout selection or fetch its pinned SDSS spectra."""

import argparse
import csv
import hashlib
import io
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import urlopen

from open_space import archive, cohort
from open_space.dataio import load_observation


TARGETS = Path(__file__).resolve().parents[1] / "data" / "holdout-targets.csv"
DISCOVERY_VISIT = (6497, 56329)
ID_COLUMNS = ("plate", "mjd", "fiber", "ra", "dec")
PLATE_SQL = """SELECT TOP 256 plate,mjd,ra,dec,run2d,programname,platequality,survey FROM PlateX WHERE plate BETWEEN 6497 AND 6752 AND survey='boss' AND programname='boss' AND platequality='good' AND run2d='v5_13_2' ORDER BY plate,mjd"""
TARGET_SQL = """SELECT s.plate,s.mjd,s.fiberID AS fiber,s.ra,s.dec
FROM SpecObjAll AS s
WHERE s.plate BETWEEN 6497 AND 6560 AND s.sciencePrimary=1
  AND s.survey='boss' AND s.programname='boss'
  AND s.class='GALAXY' AND s.class_noqso='GALAXY'
  AND s.zWarning_noqso=0 AND s.z_noqso BETWEEN 0.02 AND 0.10
  AND s.run2d='v5_13_2'
ORDER BY s.plate,s.mjd,s.fiberID
"""


def coordinates(row):
    ra, dec = float(row["ra"]), float(row["dec"])
    if not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError("Invalid celestial coordinates")
    return ra, dec


def distance(a, b):
    """Great-circle separation in degrees, including coordinates across RA zero."""
    ra1, dec1 = map(math.radians, a)
    ra2, dec2 = map(math.radians, b)
    half_chord_squared = (math.sin((dec1-dec2)/2)**2 +
                          math.cos(dec1)*math.cos(dec2)*math.sin((ra1-ra2)/2)**2)
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, half_chord_squared))))


def visit(row):
    return int(row["plate"]), int(row["mjd"])


def select_fields(plate_rows, target_rows):
    """Freeze three separated visits and up to eight distinct targets per visit."""
    plate_rows = sorted(plate_rows, key=visit)
    galaxies = {}
    for row in target_rows:
        galaxies.setdefault(visit(row), []).append(row)
    discovery = next((r for r in plate_rows if visit(r) == DISCOVERY_VISIT), None)
    if not discovery or not galaxies.get(DISCOVERY_VISIT):
        raise ValueError("Missing discovery center or catalog target positions")
    centers = [coordinates(discovery)]
    prior_targets = [coordinates(r) for r in galaxies[DISCOVERY_VISIT]]
    used_dates = {DISCOVERY_VISIT[1]}
    used_plates = {DISCOVERY_VISIT[0]}
    selected = []
    for field in plate_rows:
        plate, mjd = visit(field)
        if plate <= DISCOVERY_VISIT[0] or plate > 6560 or plate in used_plates or mjd in used_dates:
            continue
        center = coordinates(field)
        if any(distance(center, other) < 8.0 for other in centers):
            continue
        group = []
        positions = []
        for target in sorted(galaxies.get((plate, mjd), ()), key=lambda r: int(r["fiber"])):
            pos = coordinates(target)
            if any(distance(pos, other) <= 2.0/3600 for other in prior_targets + positions):
                continue
            group.append(target)
            positions.append(pos)
            if len(group) == 8:
                break
        if len(group) < 6:
            continue
        selected.extend(group)
        prior_targets.extend(positions)
        centers.append(center)
        used_dates.add(mjd)
        used_plates.add(plate)
        if len(used_plates) == 4:
            return selected
    raise ValueError("Fewer than three eligible independent fields in the catalog window")


def read_catalog(sql):
    url = cohort.CATALOG + "?" + urlencode({"cmd": sql, "format": "csv"})
    with urlopen(url, timeout=60) as response:
        raw = cohort.read_bounded(response)
    records = csv.DictReader(line for line in raw.decode("utf-8-sig").splitlines()
                             if not line.startswith("#"))
    return raw, list(records)


def selection_csv(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=ID_COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row[key] for key in ID_COLUMNS})
    return stream.getvalue().encode()


def science_exposures(raw):
    """Return unique exposure numbers, collapsing cameras for the same integration."""
    from astropy.io import fits

    with fits.open(io.BytesIO(raw)) as hdus:
        header = hdus[0].header
        exposures = set()
        for index in range(1, int(header["NEXP"]) + 1):
            value = str(header[f"EXPID{index:03}"])
            match = re.fullmatch(r"[br][12]-(\d{8})-\d{8}-\d{8}", value.strip())
            if not match:
                raise ValueError(f"Unknown SDSS EXPID convention: {value}")
            exposures.add(match.group(1))
    if not exposures:
        raise ValueError("No science exposure ID in the FITS header")
    return sorted(exposures)


def fetch_pinned(rows, output, fits_dir=None):
    """Use pinned FITS bytes and check the published derivative and exposure IDs."""
    output.mkdir(parents=True, exist_ok=True)
    for row in rows:
        observation = (*visit(row), int(row["fiber"]))
        key, url = archive.archive_source(*observation)
        sidecar = output / f"sdss-{key}.json"
        cached = output / "raw" / f"{key}.fits"
        # An explicitly supplied source directory allows full offline reconstruction.
        local = (fits_dir / cached.name) if fits_dir is not None else (cached if cached.is_file() else None)
        raw, source_hash = archive.read_source(url, local, row["source_sha256"])
        if " ".join(science_exposures(raw)) != row["science_exposures"]:
            raise ValueError(f"Science exposure mismatch for {key}")
        cached.parent.mkdir(exist_ok=True)
        if local != cached:
            archive.write_immutable(cached, raw)
        if not sidecar.exists() or not sidecar.with_suffix(".csv").exists():
            archive.ingest(observation, output, local_fits=cached, expected_hash=source_hash)
        if hashlib.sha256(sidecar.read_bytes()).hexdigest() != row["manifest_sha256"]:
            raise ValueError(f"Observation manifest mismatch for {key}")
        manifest, _ = load_observation(sidecar)
        if (manifest["observation_id"] != f"sdss-{key}" or
                any(manifest[k] != row[k] for k in ("source_sha256", "csv_sha256"))):
            raise ValueError(f"Provenance mismatch for {key}")
        if (manifest["plate_quality"] != "good" or manifest["program"] != "boss" or
                manifest["science_primary"] != 1 or manifest["classification"] != "GALAXY" or
                manifest["redshift_column"] != "SPALL.Z_NOQSO" or manifest["zwarning"] != 0 or
                not 0.02 <= manifest["redshift"] <= 0.10):
            raise ValueError(f"Observation fails FITS eligibility checks: {key}")
        print(key, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    tasks = parser.add_subparsers(dest="task", required=True)
    check = tasks.add_parser("select", help="Requery DR17 metadata and compare selected targets")
    check.add_argument("--output", required=True, type=Path)
    fetch = tasks.add_parser("fetch", help="Download or verify the pinned observing sample")
    fetch.add_argument("--output", required=True, type=Path)
    fetch.add_argument("--fits-dir", type=Path, help="Offline directory of PLATE-MJD-FIBER.fits files")
    args = parser.parse_args()
    with TARGETS.open() as stream:
        pinned = list(csv.DictReader(stream))
    if args.task == "select":
        plate_raw, plates = read_catalog(PLATE_SQL)
        galaxy_raw, galaxies = read_catalog(TARGET_SQL)
        actual = selection_csv(select_fields(plates, galaxies))
        if actual != selection_csv(pinned):
            raise ValueError("Live metadata selection differs from frozen target list")
        args.output.mkdir(parents=True, exist_ok=True)
        for name, raw in (("plates.csv", plate_raw), ("galaxies.csv", galaxy_raw),
                          ("frozen-ids.csv", actual)):
            archive.write_immutable(args.output / name, raw)
            print(name, hashlib.sha256(raw).hexdigest())
    else:
        fetch_pinned(pinned, args.output, args.fits_dir)


if __name__ == "__main__":
    main()
