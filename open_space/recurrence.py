"""Rank recurrent unidentified features across distinct sky positions in evidence packets."""

import argparse
from bisect import bisect_right
from collections import defaultdict
import json
import math
from pathlib import Path
from statistics import median

from .dataio import canonical, digest


C_KMS = 299792.458
ARCSEC = math.pi / (180 * 3600)


def packet_files(inputs):
    """Accept explicit reports, a report directory or a campaign containing reports/."""
    found = set()
    for path in inputs:
        if path.is_file():
            found.add(path.resolve())
        elif path.is_dir():
            nested = sorted(path.glob("**/reports/*.json"))
            found.update(p.resolve() for p in (nested or sorted(path.glob("*.json"))))
        else:
            raise ValueError(f"Missing input: {path}")
    if not found:
        raise ValueError("No evidence packets found")
    return sorted(found)


def _target_groups(observations, radius):
    """Join positional neighbors conservatively, including RA wrap and both poles."""
    # Chord distance is monotonic with angular separation. A 3D grid avoids an
    # all-pairs search. Connected components merge close chains conservatively.
    side = 2 * math.sin(radius * ARCSEC / 2)
    cells = defaultdict(list)
    parents = list(range(len(observations)))

    def root(i):
        while i != parents[i]:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i

    for i, obs in enumerate(observations):
        ra, dec = float(obs["ra_deg"]), float(obs["dec_deg"])
        if not 0 <= ra < 360 or not -90 <= dec <= 90:
            raise ValueError("Invalid sky coordinates")
        ra, dec = math.radians(ra), math.radians(dec)
        xyz = (math.cos(ra) * math.cos(dec), math.sin(ra) * math.cos(dec), math.sin(dec))
        x, y, z = cell = tuple(math.floor(p / side) for p in xyz)
        for a in (x - 1, x, x + 1):
            for b in (y - 1, y, y + 1):
                for c in (z - 1, z, z + 1):
                    for j, other in cells.get((a, b, c), ()):
                        if math.dist(xyz, other) <= side:
                            parents[root(i)] = root(j)
        cells[cell].append((i, xyz))
    return [root(i) for i in range(len(observations))]


def _clusters(hits, wave_key, tolerance, minimum):
    result = []
    factor = math.exp(tolerance / C_KMS)
    for direction in ("absorption", "emission"):
        ordered = sorted((h for h in hits if h["direction"] == direction),
                         key=lambda h: (h[wave_key], h["observation_id"], h["feature_id"]))
        waves = [h[wave_key] for h in ordered]
        windows = []
        last_end = 0
        for start, wave in enumerate(waves):
            end = bisect_right(waves, wave * factor)
            # Same right endpoint makes this a subset of the preceding window.
            if end == last_end:
                continue
            last_end = end
            members = {}
            for item in ordered[start:end]:
                other = members.get(item["target"])
                if other is None or item["score"] > other["score"]:
                    members[item["target"]] = item
            if len(members) >= minimum:
                reps = sorted(members.values(), key=lambda h: (h[wave_key], h["observation_id"], h["feature_id"]))
                width = C_KMS * math.log(reps[-1][wave_key] / reps[0][wave_key])
                strength = median(h["score"] for h in reps)
                windows.append((-len(reps), -strength, width, start, end, reps))
        # Overlapping windows describe competing memberships. Prefer the one with
        # more distinct targets, then stronger median rank and then smaller span.
        taken = set()
        for neg_count, neg_strength, width, start, end, reps in sorted(windows):
            if any(i in taken for i in range(start, end)):
                continue
            taken.update(range(start, end))
            observed = [h["wavelength"] for h in reps]
            rest = [h["rest_wavelength"] for h in reps]
            visits = sorted({(h["plate"], h["mjd"]) for h in reps})
            result.append({
                "direction": direction, "target_count": -neg_count, "visit_count": len(visits),
                "visits": [{"plate": p, "mjd": m} for p, m in visits],
                "median_score": -neg_strength, "frame_span_kms": width,
                "median_rest_angstrom": median(rest), "median_observed_angstrom": median(observed),
                "observed_span_angstrom": max(observed) - min(observed),
                "rest_span_angstrom": max(rest) - min(rest),
                "members": [{k: v for k, v in h.items() if k != "target"} for h in reps],
            })
    return sorted(result, key=lambda g: (-g["target_count"], -g["median_score"],
                                        g["frame_span_kms"], g["median_rest_angstrom"], g["direction"]))


def summarize(packets, tolerance_kms=150.0, min_targets=2, target_radius_arcsec=2.0):
    if not 0 < tolerance_kms <= 10000 or not 0 < target_radius_arcsec <= 60 or min_targets < 2:
        raise ValueError("Require 0 < tolerance <= 10000 km/s, 0 < radius <= 60 arcsec and at least 2 targets")
    compatibility = None
    seen = {}
    identifiers = {}
    skipped, duplicates = 0, 0
    for packet in packets:
        if packet["schema_version"] != 1 or packet["claim_status"] != "screening_only":
            raise ValueError("Expected an Open-Space evidence packet")
        if packet["state"] != "completed":
            skipped += 1
            continue
        scan, obs = packet["stages"]["detection"], packet["observation"]
        if (scan["source_context"] != {"kind": "observed", "observation_id": obs["observation_id"]}
                or obs["wavelength_frame"] != "observed vacuum" or obs["wavelength_unit"] != "Angstrom"):
            raise ValueError("Expected a traceable observation on a vacuum Angstrom grid")
        settings = {"software_sha256": packet["software_sha256"], "detector_schema": scan["schema_version"],
                    "method": scan["method"], "catalog": scan["catalog"], "score_threshold": scan["score_threshold"]}
        if compatibility is not None and settings != compatibility:
            raise ValueError("Incompatible scans: separate evidence by code, detector, catalog and threshold")
        compatibility = settings
        source = obs["source_sha256"]
        key = (obs["observation_id"], source)
        # Re-running or changing a review cannot create observational evidence.
        # Conflicting copies of a coadd must be reconciled upstream.
        signature = digest(canonical([scan, obs]).encode())
        if source in seen:
            if seen[source][0] != signature:
                raise ValueError(f"Conflicting evidence for source {source}")
            duplicates += 1
            continue
        if key[0] in identifiers:
            raise ValueError(f"Conflicting source bytes for {key[0]}")
        identifiers[key[0]] = source
        seen[source] = (signature, obs, scan)

    records = sorted(seen.values(), key=lambda rec: rec[1]["observation_id"])
    groups = _target_groups([r[1] for r in records], target_radius_arcsec)
    hits = []
    for target, (_, obs, scan) in zip(groups, records):
        for i, feature in enumerate(scan["features"], 1):
            if feature["classification"] != "unidentified":
                continue
            h = {k: feature[k] for k in ("direction", "wavelength", "rest_wavelength", "score")}
            if (h["direction"] not in ("absorption", "emission") or
                    not all(math.isfinite(h[k]) and h[k] > 0 for k in ("wavelength", "rest_wavelength", "score"))):
                raise ValueError("Invalid unidentified feature")
            h.update(target=target, feature_id=f"F{i:04d}", observation_id=obs["observation_id"],
                     plate=obs["plate"], mjd=obs["mjd"], source_sha256=obs["source_sha256"])
            hits.append(h)
    return {
        "schema_version": 1, "claim_status": "screening_only", "compatibility": compatibility,
        "tolerance_kms": tolerance_kms, "min_targets": min_targets, "target_radius_arcsec": target_radius_arcsec,
        "unique_observations": len(records), "distinct_targets": len(set(groups)),
        "skipped_incomplete": skipped, "duplicate_completed_packets": duplicates,
        "unidentified_features": len(hits),
        "rest_frame_groups": _clusters(hits, "rest_wavelength", tolerance_kms, min_targets),
        "observer_frame_control_groups": _clusters(hits, "wavelength", tolerance_kms, min_targets),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Evidence files, reports directories or campaign roots")
    parser.add_argument("--tolerance-kms", type=float, default=150)
    parser.add_argument("--min-targets", type=int, default=2)
    parser.add_argument("--target-radius-arcsec", type=float, default=2)
    args = parser.parse_args(argv)
    try:
        files = packet_files(args.inputs)
        result = summarize((json.loads(p.read_bytes()) for p in files),
                           args.tolerance_kms, args.min_targets, args.target_radius_arcsec)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (OSError, KeyError, TypeError, ValueError) as exc:
        parser.exit(1, f"recurrence failed: {exc}\n")


if __name__ == "__main__":
    main()
