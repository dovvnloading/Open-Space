"""Fail CI on missing axiom reports or dependencies outside standard Lean."""
import re
import sys
from pathlib import Path

GENERAL = {
    "commonShiftPreservesDuration", "commonShiftPreservesWindow", "orderedOverlap",
    "orderedCoveragePartition", "largerLimbPreservesCoverage", "calibrationAccounting",
    "noAssessableIff", "pairedOutcomePartition", "stellarRadiusScaling",
}
EXACT = {
    "nicerClockOffset", "nicer2018PrintedLabels", "nicer2022PrintedLabels",
    "nicerCatalogExact", "decemberContactLabels", "decemberContactOrder",
    "decemberContactDurations", "decemberDurationRounding",
    *(f"oecAxis{i}" for i in range(6)),
}
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}
PATTERN = re.compile(
    r"^'OpenSpace\.([A-Za-z0-9_]+)' "
    r"(?:does not depend on any axioms|depends on axioms: \[([^\]]*)\])$"
)


def main(paths):
    seen = {}
    for path in paths:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            match = PATTERN.fullmatch(line.strip())
            if not match:
                continue
            name, dependency_list = match.groups()
            if name in seen:
                raise ValueError(f"Duplicate axiom report: {name}")
            seen[name] = {item.strip() for item in (dependency_list or "").split(",")
                          if item.strip()}
    expected = GENERAL | EXACT
    if set(seen) != expected:
        raise ValueError(f"Missing/unexpected axiom reports: {expected - set(seen)}, {set(seen) - expected}")
    for name, dependencies in seen.items():
        allowed = set() if name in EXACT else ALLOWED
        if dependencies - allowed:
            raise ValueError(f"{name} uses forbidden dependencies: {dependencies - allowed}")
    print(f"Axiom audit passed for {len(seen)} theorems; {len(EXACT)} exact certificates have no axioms.")


if __name__ == "__main__":
    main(sys.argv[1:])
