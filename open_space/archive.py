"""Fetch or replay an SDSS DR17 BOSS/eBOSS lite spectrum."""

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
from urllib.error import URLError
from urllib.request import urlopen


PILOT = (6497, 56329, 139)
PILOT_SHA256 = "5a06bc1f02dd76b17323e5888e22149e315b10b3a2f60c5a47dbffa19d391704"
RUN2D = "v5_13_2"
ARCHIVE = "https://data.sdss.org/sas/dr17/eboss/spectro/redux"
MAX_BYTES = 8 * 1024 * 1024


def archive_source(plate, mjd, fiber):
    if not (1 <= plate <= 99999 and 10000 <= mjd <= 99999 and 1 <= fiber <= 1000):
        raise ValueError("Invalid DR17 plate, five-digit MJD, or fiber (1 to 1000)")
    key = f"{plate:04d}-{mjd:05d}-{fiber:04d}"
    return key, f"{ARCHIVE}/{RUN2D}/spectra/lite/{plate:04d}/spec-{key}.fits"


def read_source(source_url, local_fits, expected_hash):
    if expected_hash and not re.fullmatch(r"[a-fA-F0-9]{64}", expected_hash):
        raise ValueError("Expected SHA-256 must contain 64 hexadecimal characters")
    stream = local_fits.open("rb") if local_fits else urlopen(source_url, timeout=30)
    with stream:
        raw = stream.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise ValueError("Spectrum exceeds the supported DR17 lite product size")
    digest = hashlib.sha256(raw).hexdigest()
    if expected_hash and digest != expected_hash.lower():
        raise ValueError(f"Source checksum mismatch. Retrieved SHA-256 {digest}")
    return raw, digest


def csv_and_metadata(raw, observation):
    try:
        from astropy.io import fits
    except ModuleNotFoundError as exc:
        raise ValueError("Install Astropy to convert FITS files: python -m pip install astropy") from exc

    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("wavelength", "flux", "ivar", "mask", "and_mask", "or_mask"))
    usable = 0
    previous = 0.0
    with fits.open(io.BytesIO(raw)) as hdus:
        header = hdus[0].header
        if header.get("VACUUM") is not True or str(header.get("RUN2D", "")).strip() != RUN2D:
            raise ValueError(f"Expected SDSS vacuum wavelengths and pipeline {RUN2D}")
        if tuple(header.get(k) for k in ("PLATEID", "MJD", "FIBERID")) != observation:
            raise ValueError("FITS primary header does not match requested observation")
        if len(hdus["SPALL"].data) != 1:
            raise ValueError("Expected one pipeline fit for this observation")
        obj = hdus["SPALL"].data[0]
        if tuple(int(obj[k]) for k in ("PLATE", "MJD", "FIBERID")) != observation:
            raise ValueError("FITS pipeline fit does not match requested observation")
        classification = str(obj["CLASS"]).strip()
        suffix = "_NOQSO" if classification == "GALAXY" else ""
        redshift = float(obj[f"Z{suffix}"])
        if not math.isfinite(redshift) or redshift <= -1:
            raise ValueError("Source has no usable redshift")
        coadd = hdus["COADD"].data
        for loglam, flux, ivar, and_mask, or_mask in zip(
            coadd["loglam"], coadd["flux"], coadd["ivar"], coadd["and_mask"], coadd["or_mask"], strict=True
        ):
            wave = 10 ** float(loglam)
            if not math.isfinite(wave) or wave <= previous:
                raise ValueError("Wavelengths must be finite, positive, and strictly increasing")
            previous = wave
            flux, ivar = float(flux), float(ivar)
            if not math.isfinite(flux):
                flux, ivar = 0.0, 0.0
            if not math.isfinite(ivar) or ivar <= 0:
                ivar = 0.0
            # SDSS stores these flags in a signed 32-bit FITS integer.
            always_flagged = int(and_mask) & 0xffffffff
            sometimes_flagged = int(or_mask) & 0xffffffff
            usable += ivar > 0 and always_flagged == 0
            writer.writerow((f"{wave:.8f}", f"{flux:.9g}", f"{ivar:.9g}", always_flagged, always_flagged, sometimes_flagged))
        if not usable:
            raise ValueError("Spectrum has no usable pixels after pixel-quality filtering")
        redshift_error = float(obj[f"Z_ERR{suffix}"])
        metadata = {
            "redshift": redshift,
            "redshift_error": redshift_error if math.isfinite(redshift_error) and redshift_error >= 0 else None,
            "redshift_column": f"SPALL.Z{suffix}",
            "zwarning": int(obj[f"ZWARNING{suffix}"]),
            "classification": str(obj[f"CLASS{suffix}"]).strip(),
            "subclass": str(obj[f"SUBCLASS{suffix}"]).strip(),
            "ra_deg": float(obj["PLUG_RA"]),
            "dec_deg": float(obj["PLUG_DEC"]),
            "plate_quality": str(obj["PLATEQUALITY"]).strip(),
            "program": str(obj["PROGRAMNAME"]).strip(),
            "science_primary": int(obj["SPECPRIMARY"]),
            "pixels": len(coadd),
            "usable_pixels": usable,
        }
    return output.getvalue().encode("ascii"), metadata


def write_immutable(path, content):
    """Publish a completed file atomically; an identical replay is harmless."""
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as staged:
        staged_path = Path(staged.name)
    try:
        staged_path.write_bytes(content)
        try:
            os.link(staged_path, path)
        except FileExistsError:
            if path.read_bytes() != content:
                raise ValueError(f"Conflicting output {path}. Choose a new output directory")
    finally:
        staged_path.unlink(missing_ok=True)


def ingest(observation, output_dir, local_fits=None, expected_hash=None):
    key, url = archive_source(*observation)
    if observation == PILOT:
        if expected_hash and expected_hash.lower() != PILOT_SHA256:
            raise ValueError("Expected SHA-256 conflicts with the pinned pilot source")
        expected_hash = PILOT_SHA256
    raw, raw_hash = read_source(url, local_fits, expected_hash)
    csv_bytes, metadata = csv_and_metadata(raw, observation)
    stem = f"sdss-{key}"
    manifest = {
        "schema_version": 1,
        "observation_id": stem,
        "archive": "Sloan Digital Sky Survey",
        "release": "DR17",
        "product": "BOSS/eBOSS lite coadd",
        "run2d": RUN2D,
        "plate": observation[0],
        "mjd": observation[1],
        "fiber": observation[2],
        "source_url": url,
        "source_sha256": raw_hash,
        "source_bytes": len(raw),
        "csv_file": f"{stem}.csv",
        "csv_sha256": hashlib.sha256(csv_bytes).hexdigest(),
        "wavelength_unit": "Angstrom",
        "wavelength_frame": "observed vacuum",
        "flux_unit": "1e-17 erg cm^-2 s^-1 Angstrom^-1",
        "mask_column": "COADD.and_mask",
        "eligibility": "ivar > 0 and mask == 0; or_mask retained for follow-up",
        **metadata,
    }
    json_bytes = (json.dumps(manifest, indent=2, allow_nan=False) + "\n").encode()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_immutable(output_dir / manifest["csv_file"], csv_bytes)
    manifest_path = output_dir / f"{stem}.json"
    write_immutable(manifest_path, json_bytes)
    return manifest_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="Destination directory")
    parser.add_argument("--fits", type=Path, help="Read the selected observation's lite FITS locally")
    parser.add_argument("--plate", type=int)
    parser.add_argument("--mjd", type=int)
    parser.add_argument("--fiber", type=int)
    parser.add_argument("--sha256", help="Verify source bytes against this SHA-256")
    args = parser.parse_args()
    values = (args.plate, args.mjd, args.fiber)
    if any(v is not None for v in values) and any(v is None for v in values):
        parser.error("Supply --plate, --mjd, and --fiber together")
    observation = PILOT if args.plate is None else values
    try:
        path = ingest(observation, args.output, args.fits, args.sha256)
    except (OSError, ValueError, KeyError, URLError) as exc:
        parser.exit(1, f"Archive ingestion failed: {exc}\n")
    print(path)


if __name__ == "__main__":
    main()
