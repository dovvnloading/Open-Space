"""Plan a bounded DR17 galaxy cohort and fetch its observations."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import urlopen

from .archive import PILOT, RUN2D, archive_source, ingest, write_immutable


CATALOG = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
MAX_TARGETS = 100
MAX_PLAN_BYTES = 256 * 1024


def selection_query(plate, mjd, after_fiber, limit):
    archive_source(plate, mjd, 1)
    if not 0 <= after_fiber <= 1000 or not 1 <= limit <= MAX_TARGETS:
        raise ValueError("Use a fiber cursor from 0 to 1000 and a limit from 1 to 100")
    return f"""SELECT TOP {limit} s.plate, s.mjd, s.fiberID AS fiber
FROM SpecObjAll AS s
WHERE s.plate = {plate} AND s.mjd = {mjd} AND s.fiberID > {after_fiber}
  AND s.sciencePrimary = 1 AND s.survey = 'boss' AND s.programname = 'boss'
  AND s.class = 'GALAXY' AND s.class_noqso = 'GALAXY'
  AND s.zWarning_noqso = 0 AND s.z_noqso BETWEEN 0.02 AND 0.10
  AND s.run2d = '{RUN2D}'
ORDER BY s.fiberID"""


def target_list(raw, plate, mjd, after_fiber, limit):
    lines = (line for line in raw.decode("utf-8-sig").splitlines() if not line.startswith("#"))
    records = csv.DictReader(lines)
    if not {"plate", "mjd", "fiber"}.issubset(records.fieldnames or ()):
        raise ValueError("Expected CSV header plate,mjd,fiber from the DR17 catalog or operator")
    targets = []
    seen = set()
    for row in records:
        observation = tuple(int(row[k]) for k in ("plate", "mjd", "fiber"))
        archive_source(*observation)
        if observation[:2] != (plate, mjd) or observation[2] <= after_fiber:
            raise ValueError("Target is outside this plate, MJD, or fiber cursor")
        if observation in seen:
            raise ValueError("Duplicate observation in cohort")
        seen.add(observation)
        target = dict(zip(("plate", "mjd", "fiber"), observation))
        digest = (row.get("source_sha256") or "").strip().lower()
        if digest:
            if not re.fullmatch(r"[a-f0-9]{64}", digest):
                raise ValueError("Invalid source_sha256 in target list")
            target["source_sha256"] = digest
        targets.append(target)
        if len(targets) > limit:
            raise ValueError("Target list exceeds the requested limit")
    return sorted(targets, key=lambda t: t["fiber"])


def encode(document):
    return (json.dumps(document, indent=2, allow_nan=False) + "\n").encode()


def read_bounded(stream):
    raw = stream.read(MAX_PLAN_BYTES + 1)
    if len(raw) > MAX_PLAN_BYTES:
        raise ValueError("Target plan or catalog response exceeds the expected size")
    return raw


def plan(plate, mjd, after_fiber, limit, input_csv, output):
    query = selection_query(plate, mjd, after_fiber, limit)
    if input_csv:
        origin = "operator_csv"
        with input_csv.open("rb") as stream:
            raw = read_bounded(stream)
    else:
        origin = "skyserver"
        url = CATALOG + "?" + urlencode({"cmd": query, "format": "csv"})
        with urlopen(url, timeout=60) as stream:
            raw = read_bounded(stream)
    targets = target_list(raw, plate, mjd, after_fiber, limit)
    doc = {
        "schema_version": 1, "release": "DR17", "origin": origin,
        "catalog_url": CATALOG if not input_csv else None,
        "query": query if not input_csv else None,
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "plate": plate, "mjd": mjd, "after_fiber": after_fiber,
        "limit": limit, "targets": targets,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    write_immutable(output, encode(doc))
    return doc


def fetch(plan_path, output):
    with plan_path.open("rb") as stream:
        plan_bytes = read_bounded(stream)
    doc = json.loads(plan_bytes)
    if doc["schema_version"] != 1 or doc["release"] != "DR17":
        raise ValueError("Unsupported cohort plan schema or data release")
    plate, mjd, cursor, limit = (doc[k] for k in ("plate", "mjd", "after_fiber", "limit"))
    selection_query(plate, mjd, cursor, limit)
    targets = doc["targets"]
    if not isinstance(targets, list) or len(targets) > limit:
        raise ValueError("Plan exceeds requested cohort size")
    previous = cursor
    # Check every identifier before making network requests.
    for target in targets:
        p, m, f = (target[k] for k in ("plate", "mjd", "fiber"))
        archive_source(p, m, f)
        if (p, m) != (plate, mjd) or f <= previous:
            raise ValueError("Plan must contain ordered unique targets from one visit")
        previous = f
    output.mkdir(parents=True, exist_ok=True)
    cohort = []
    for target in targets:
        obs = (plate, mjd, target["fiber"])
        key, _ = archive_source(*obs)
        expected = target.get("source_sha256")
        # Existing manifests pin the raw digest before re-fetch. A rerun compares
        # both products byte-for-byte; partial runs can finish after interruptions.
        previous_path = output / f"sdss-{key}.json"
        if previous_path.is_file():
            seen_hash = json.loads(previous_path.read_bytes())["source_sha256"]
            if expected and seen_hash != expected:
                raise ValueError(f"Previously recorded source differs for {key}")
            expected = seen_hash
        sidecar = ingest(obs, output, expected_hash=expected)
        metadata = json.loads(sidecar.read_bytes())
        if (metadata["program"] != "boss" or metadata["science_primary"] != 1
                or metadata["redshift_column"] != "SPALL.Z_NOQSO"
                or metadata["classification"] != "GALAXY" or metadata["zwarning"] != 0
                or not 0.02 <= metadata["redshift"] <= 0.10):
            raise ValueError(f"{key} fails the galaxy cohort criteria in its FITS metadata")
        cohort.append({
            "observation_id": metadata["observation_id"], "provenance": sidecar.name,
            "csv_file": metadata["csv_file"], "redshift": metadata["redshift"],
            "source_sha256": metadata["source_sha256"], "csv_sha256": metadata["csv_sha256"],
        })
    digest = hashlib.sha256(plan_bytes).hexdigest()
    receipt = output / f"cohort-{digest[:16]}.json"
    write_immutable(receipt, encode({"plan_sha256": digest, "observations": cohort}))
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    tasks = parser.add_subparsers(dest="command", required=True)
    prepare = tasks.add_parser("plan", help="Freeze IDs from one plate visit")
    prepare.add_argument("--plate", type=int, default=PILOT[0])
    prepare.add_argument("--mjd", type=int, default=PILOT[1])
    prepare.add_argument("--after-fiber", type=int, default=0)
    prepare.add_argument("--limit", type=int, default=10)
    prepare.add_argument("--input-csv", type=Path)
    prepare.add_argument("--output", type=Path, required=True)
    download = tasks.add_parser("fetch", help="Verify and normalize a saved cohort")
    download.add_argument("--plan", type=Path, required=True)
    download.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "plan":
            doc = plan(args.plate, args.mjd, args.after_fiber, args.limit, args.input_csv, args.output)
            print(f'{len(doc["targets"])} targets in {args.output}')
        else:
            print(fetch(args.plan, args.output))
    except (ValueError, OSError, KeyError, TypeError) as exc:
        parser.exit(1, f"Cohort {args.command} failed: {exc}\n")


if __name__ == "__main__":
    main()
