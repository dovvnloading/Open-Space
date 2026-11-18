"""Verify PDS archives and extract complete paired Alice histogram windows."""

from argparse import ArgumentParser
from collections import Counter
from hashlib import file_digest, sha256
import io
import json
from pathlib import Path
import tarfile

from astropy.io import fits
import numpy as np


FIELDS = (
    tuple(f"NCOB{number:02d}" for number in range(1, 16))
    + ("DCAL01", "DCAL02", "DCAL03", "DCAL05", "DCAL06", "DCAL07", "DCAL08", "SHOCK1", "SHOCK2")
)
NORTH = {"NCOB12", "NCOB13", "NCOB15", "DCAL03", "DCAL05", "DCAL07", "DCAL08", "SHOCK1", "SHOCK2"}
CENTERS = (830.0, 1134.6, 1549.5)
SUPPORT = 60.0


def read_archive(path, source, calibrated):
    with path.open("rb") as stream:
        digest = file_digest(stream, "sha256").hexdigest()
    if path.stat().st_size != source["bytes"] or digest != source["raw_sha256"]:
        raise ValueError(f"PDS archive digest or size mismatch: {path.name}")
    selected = {}
    rejected = []
    with tarfile.open(path, "r|gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".fit"):
                continue
            content = archive.extractfile(member).read()
            with fits.open(io.BytesIO(content), memmap=False) as hdus:
                header = hdus[0].header
                field = header.get("SAPDESC", "").removeprefix("ALICE CUVB ")
                if field not in FIELDS:
                    continue
                duration = float(header["EXPTIME"])
                door = int(header["APRDRPOS"])
                start = float(header["STARTMET"])
                stop = float(header["STOPMET"])
                mode = int(header["ACQMODE"])
                key = (field, start, stop, door)
                if (duration != 3600 or mode != 1 or door not in (1, 2)
                        or header["OBSCOMPL"] != "COMPLETE" or stop - start != duration):
                    rejected.append({"path": member.name, "reason": "incompatible histogram exposure"})
                    continue
                if key in selected:
                    raise ValueError(f"Duplicate Alice observation key: {key}")
                array = hdus["WAVELENGTH_IMAGE"].data if calibrated else hdus[0].data
                if array.shape != (32, 1024):
                    raise ValueError(f"Unexpected Alice histogram shape: {member.name}")
                if not calibrated and (array.dtype.kind not in "iu" or np.min(array) < 0):
                    raise ValueError(f"Invalid raw histogram counts: {member.name}")
                record = {
                    "path": member.name, "sha256": sha256(content).hexdigest(),
                    "met": int(header["MET"]), "visit": int(header["VISITID"]),
                    "field": field, "start": start, "stop": stop, "door": door,
                    "duration": duration, "hv_safety": int(header["HVSFTY"]),
                    "hv_mask": int(header["HVSMSK"]),
                }
                selected[key] = record, array.astype(np.float32 if calibrated else np.int32)
    return selected, rejected


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("sources", type=Path)
    parser.add_argument("archive_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.sources.read_bytes())
    raw, raw_rejected = read_archive(args.archive_dir / manifest["raw"]["file"], manifest["raw"], False)
    cal, cal_rejected = read_archive(args.archive_dir / manifest["calibrated"]["file"], manifest["calibrated"], True)
    by_field = Counter(key[0] for key in raw.keys() & cal.keys())
    if by_field.keys() != set(FIELDS):
        raise ValueError(f"Missing fields: {set(FIELDS) - by_field.keys()}")
    columns = np.zeros(1024, dtype=bool)
    paired = sorted(raw.keys() & cal.keys())
    for key in paired:
        wave = cal[key][1]
        near = (wave >= 1195) & (wave <= 1235)
        for center in CENTERS:
            near |= np.abs(wave - center) <= SUPPORT
        columns |= near.any(axis=0)
    records = []
    histograms, wavelengths = [], []
    for key in paired:
        raw_meta, histogram = raw[key]
        cal_meta, wave = cal[key]
        if any(raw_meta[field] != cal_meta[field] for field in ("met", "visit", "hv_safety", "hv_mask")):
            raise ValueError(f"Raw/calibrated metadata mismatch: {key}")
        record = dict(raw_meta, calibrated_path=cal_meta["path"], calibrated_sha256=cal_meta["sha256"])
        record["hemisphere"] = "north" if key[0] in NORTH else "south"
        records.append(record)
        histograms.append(histogram[:, columns])
        wavelengths.append(wave[:, columns])
    assert len(records) == len(histograms) == len(wavelengths)
    args.output_dir.mkdir()
    np.savez_compressed(args.output_dir / "windows.npz", counts=histograms, wave=wavelengths,
                       columns=np.flatnonzero(columns))
    output = {
        "records": records, "column_count": int(columns.sum()), "fields": dict(sorted(by_field.items())),
        "raw_only": [raw[key][0] for key in sorted(raw.keys() - cal.keys())],
        "calibrated_only": [cal[key][0] for key in sorted(cal.keys() - raw.keys())],
        "rejected": raw_rejected + cal_rejected,
    }
    with (args.output_dir / "observations.json").open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
