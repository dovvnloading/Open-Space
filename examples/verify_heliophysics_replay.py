"""Fail archival CI if source bytes or committed science results drift.

The research replays are deterministic and use pinned source digests. A changed
upstream source must be re-reviewed and accompanied by a reviewed new reference
packet, rather than silently producing a different statistic on green CI.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

SOURCE_SHA256 = {
    "forbush": {
        "feid": "6d88f455a97c5e2baf9a2d6e63f0f6b549865907fe8a2a0a0eff24e3a5a4914f",
        "omni": "e97d99bd5ea47786d2ab39dd3f65b3f2db193c48ef4c8ba992a6c5862e13662a",
    },
    "neutron": "7fcd5a46266cae549fb2658a12525c1c581a249f4a4b19f9286294fc83ad884b",
    "cme": "8794251f244af391ab275ae6e10d5ca36c37c591461ca7dd8d3d4d31d2de0eb2",
}


def mismatches(expected: Any, actual: Any, path: str = "$", max_items: int = 12) -> list[str]:
    """Report first differences, including absent/extra provenance fields."""
    found: list[str] = []

    def walk(e: Any, a: Any, p: str) -> None:
        if len(found) >= max_items:
            return
        if (isinstance(e, (int, float)) and not isinstance(e, bool)
                and isinstance(a, (int, float)) and not isinstance(a, bool)):
            # JSON numbers 1 and 1.0 are semantically equivalent; the GitHub
            # connector normalizes integral floats to integers in commits.
            if e != a:
                found.append(f"{p}: {e!r} -> {a!r}")
        elif type(e) is not type(a):
            found.append(f"{p}: type changed {type(e).__name__} -> {type(a).__name__}")
        elif isinstance(e, dict):
            for key in sorted(e.keys() | a.keys()):
                if len(found) >= max_items:
                    break
                if key not in e:
                    found.append(f"{p}.{key}: unexpected key")
                elif key not in a:
                    found.append(f"{p}.{key}: missing key")
                else:
                    walk(e[key], a[key], f"{p}.{key}")
        elif isinstance(e, list):
            if len(e) != len(a):
                found.append(f"{p}: list length {len(e)} -> {len(a)}")
            for index, (x, y) in enumerate(zip(e, a)):
                if len(found) >= max_items:
                    break
                walk(x, y, f"{p}[{index}]")
        elif e != a:
            found.append(f"{p}: {e!r} -> {a!r}")

    walk(expected, actual, path)
    return found


def verify(kind: str, reference: dict, observed: dict) -> None:
    """Check fixed source SHA-256 plus all stored output statistics."""
    if kind == "forbush":
        digest = {
            key: observed.get("source_files", {}).get(key, {}).get("sha256")
            for key in ("feid", "omni")
        }
    elif kind == "neutron":
        digest = observed.get("input_sha256")
    elif kind == "cme":
        digest = observed.get("source_sha256")
    else:
        raise ValueError(f"Unknown experiment: {kind}")

    if digest != SOURCE_SHA256[kind]:
        raise ValueError(f"Upstream source digest changed for {kind}; review required: {digest}")

    differences = mismatches(reference, observed)
    if differences:
        raise ValueError(
            "Measurement packet differs from committed reference; review required:\n"
            + "\n".join(differences)
        )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kind", choices=sorted(SOURCE_SHA256), required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--observed", type=Path, required=True)
    a = p.parse_args()
    expected = json.loads(a.reference.read_text())
    actual = json.loads(a.observed.read_text())
    verify(a.kind, expected, actual)
    print(f"PASS: {a.kind} source digest and every reference measurement match")


if __name__ == "__main__":
    main()
