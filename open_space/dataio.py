"""Load an archived observation and verify the bytes used in every experiment."""

import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def load_observation(path):
    path = Path(path)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest["schema_version"] != 1 or manifest["wavelength_frame"] != "observed vacuum":
        raise ValueError("Expected schema 1 and observed vacuum wavelengths")
    if manifest["wavelength_unit"] != "Angstrom":
        raise ValueError("Convert the wavelength grid to Angstrom before enqueueing")
    if not isinstance(manifest["observation_id"], str) or not manifest["observation_id"].strip():
        raise ValueError("Missing observation ID")
    redshift = float(manifest["redshift"])
    if not math.isfinite(redshift) or redshift <= -1:
        raise ValueError("A measured, finite redshift greater than -1 is required")
    for key in ("source_sha256", "csv_sha256"):
        if not re.fullmatch("[0-9a-f]{64}", manifest[key]):
            raise ValueError(f"Invalid {key}")
    if not manifest["source_url"]:
        raise ValueError("Missing source location")
    filename = manifest["csv_file"]
    if not filename or Path(filename).name != filename:
        raise ValueError("Place the CSV beside its observation manifest")
    raw = (path.parent / filename).read_bytes()
    if digest(raw) != manifest["csv_sha256"]:
        raise ValueError("CSV checksum mismatch")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline=""))
    if not {"wavelength", "flux", "ivar", "mask"}.issubset(reader.fieldnames or []):
        raise ValueError("Missing required CSV columns")
    rows = list(reader)
    if len(rows) != manifest["pixels"] or not rows:
        raise ValueError("Pixel count disagrees with the manifest or the spectrum is empty")
    return manifest, rows


def receipt_paths(path):
    receipt = json.loads(path.read_text())
    paths = []
    for entry in receipt["observations"]:
        filename = entry["provenance"]
        if not filename or Path(filename).name != filename:
            raise ValueError("Invalid receipt path")
        observation = path.parent / filename
        manifest = json.loads(observation.read_text())
        for key in ("observation_id", "csv_sha256", "source_sha256", "redshift", "csv_file"):
            if manifest[key] != entry[key]:
                raise ValueError(f"Receipt differs from {filename} on {key}")
        paths.append(observation)
    return paths
