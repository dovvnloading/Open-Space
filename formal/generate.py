"""Encode committed CSV inputs into Lean literals; Lean checks the claims.
No network, floating-point arithmetic, or generated proof placeholders.
The CSV parser and its transcription into Lean are outside the Lean kernel.
"""
import csv
import re
import sys
from datetime import datetime
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NICER = ROOT / "results/nicer-magnetar-microbursts/catalog-clock/windows.csv"
CONTACTS = ROOT / "results/solar-orbiter-venus-jovian-occultations/contacts.csv"
ORBITS = ROOT / "results/oec-orbital-scale-audit/primary_checks.csv"
# These are published target values, not computed from the data.
NICER_EXPECTED = [1132, 944, 843, 1, 100, 115351, 274146,
                  146, 145, 0, 1, 809, 8258]
OEC_EXPECTED = {
    "Kepler-445 b": 2143, "Kepler-445 c": 2950, "Kepler-445 d": 4194,
    "Kepler-446 b": 1585, "Kepler-446 c": 2500, "Kepler-446 d": 3527,
}


def read_rows(path, required):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not set(required) <= set(fields):
            raise ValueError(f"Missing or duplicated column: {path}")
        rows = list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError(f"Malformed CSV row: {path}")
        return rows


def nat_or_none(value):
    if value == "":
        return None
    if not re.fullmatch(r"[0-9]+", value):
        raise ValueError(f"Expected nonnegative integer, got {value!r}")
    return int(value)


def scaled(value, factor):
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", value):
        raise ValueError(f"Expected nonnegative decimal, got {value!r}")
    result = Decimal(value) * factor
    if result != result.to_integral_value():
        raise ValueError(f"Value {value!r} cannot be encoded at scale {factor}")
    return int(result)


def numeral_tuple(values):
    return "⟨" + ", ".join(map(str, values)) + "⟩"


def option(value):
    return "none" if value is None else f"some {value}"


def block_summary(rows):
    # The kernel independently checks each block against this proposed summary.
    # This is only an acceleration: final target totals are fixed above.
    result = [0] * 13
    for a, b, other in rows:
        result[0] += 1
        if a is None or b is None:
            continue
        result[1] += 1
        outcome = 2 if b > a else 3 if b == a else 4
        result[outcome] += 1
        result[5] += a
        result[6] += b
        if other:
            result[7] += 1
            result[outcome + 6] += 1
            result[11] += a
            result[12] += b
    return result


def nicer_lines():
    rows = [
        (nat_or_none(row["narrow_literal"]),
         nat_or_none(row["narrow_plus69"]),
         row["source"] != "SGR 1935+2154")
        for row in read_rows(NICER, ["source", "narrow_literal", "narrow_plus69"])
    ]
    if not rows:
        raise ValueError("Empty NICER input")
    names = []
    lines = ["-- All input rows, including those with missing coverage."]
    for index, start in enumerate(range(0, len(rows), 24)):
        block = rows[start:start + 24]
        name = f"nicerBlock{index}"
        names.append(name)
        encoded = [
            f"  ⟨{option(a)}, {option(b)}, {str(other).lower()}⟩"
            for a, b, other in block
        ]
        lines += [
            f"def {name} : List Window := [",
            ",\n".join(encoded),
            "]",
            f"theorem {name}_cert : summarize {name} = "
            f"{numeral_tuple(block_summary(block))} := by decide",
            "",
        ]
    # A binary tree of additions keeps the expression shallow.
    expressions = [f"(summarize {name})" for name in names]
    while len(expressions) > 1:
        expressions = [
            f"(addSummary {expressions[i]} {expressions[i + 1]})"
            if i + 1 < len(expressions) else expressions[i]
            for i in range(0, len(expressions), 2)
        ]
    lines += [
        f"def nicerCatalogSummary : Summary := {expressions[0]}",
        f"theorem nicerCatalogExact : nicerCatalogSummary = "
        f"{numeral_tuple(NICER_EXPECTED)} := by",
        "  unfold nicerCatalogSummary",
        "  simp only [" + ",\n    ".join(name + "_cert" for name in names)
        + "] <;> decide",
        "#print axioms nicerCatalogExact",
        "",
    ]
    return lines


def milliseconds(value, day):
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z", value):
        raise ValueError(f"Expected exact millisecond UTC label: {value}")
    stamp = datetime.strptime(value, "%Y-%m-%dT%H:%M:%S.%fZ")
    if stamp.date().isoformat() != day:
        raise ValueError(f"Contact crosses selected day: {value}")
    return (stamp.hour * 3600 + stamp.minute * 60 + stamp.second) * 1000 + stamp.microsecond // 1000


def contact_lines():
    day = "2026-12-24"
    required = ["date", "spacecraft", "direction", "venus_height_km",
                "source_radius_in_target_radii", "ingress_utc", "egress_utc",
                "duration_s"]
    selected = {}
    for row in read_rows(CONTACTS, required):
        if (row["date"] == day and row["spacecraft"] == "SolO"
                and row["direction"] in {"earth", "jupiter"}
                and row["venus_height_km"] == "0"
                and row["source_radius_in_target_radii"] == "0"):
            key = row["direction"]
            if key in selected:
                raise ValueError(f"Duplicate center contact: {key}")
            selected[key] = row
    if set(selected) != {"earth", "jupiter"}:
        raise ValueError("Missing December center contacts")
    lines = []
    for name in ("earth", "jupiter"):
        row = selected[name]
        a = milliseconds(row["ingress_utc"], day)
        b = milliseconds(row["egress_utc"], day)
        duration = scaled(row["duration_s"], 1000)
        lines += [f"def {name}StartMs : Int := {a}",
                  f"def {name}StopMs : Int := {b}",
                  f"def {name}ReportedDurationMs : Int := {duration}"]
    lines += [
        "theorem decemberContactOrder : earthStartMs ≤ jupiterStartMs ∧",
        "    jupiterStartMs ≤ earthStopMs ∧ earthStopMs ≤ jupiterStopMs := by decide",
        "theorem decemberContactDurations :",
        "    earthStopMs - earthStartMs = 139303 ∧",
        "    jupiterStopMs - jupiterStartMs = 887696 ∧",
        "    overlapMs earthStartMs earthStopMs jupiterStartMs jupiterStopMs = 18654 ∧",
        "    jupiterStartMs - earthStartMs = 120649 ∧",
        "    jupiterStopMs - earthStopMs = 869042 := by decide",
        "-- Each of two rounded endpoints contributes up to 0.5 ms;",
        "-- the rounded duration contributes up to 5 ms.",
        "theorem decemberDurationRounding :",
        "    (earthStopMs - earthStartMs - earthReportedDurationMs).natAbs ≤ 6 ∧",
        "    (jupiterStopMs - jupiterStartMs - jupiterReportedDurationMs).natAbs ≤ 6 := by decide",
        "#print axioms decemberContactDurations",
        "#print axioms decemberDurationRounding",
        "",
    ]
    return lines


def orbit_lines():
    rows = read_rows(ORBITS, ["planet", "paper_a_over_rstar",
                              "oec_host_radius_solar", "oec_reported_a_au"])
    by_name = {row["planet"]: row for row in rows}
    if len(rows) != len(by_name) or set(by_name) != set(OEC_EXPECTED):
        raise ValueError("Expected exactly six uniquely named orbital audit rows")
    lines = []
    for index, (name, rounded) in enumerate(OEC_EXPECTED.items()):
        row = by_name[name]
        q = scaled(row["paper_a_over_rstar"], 100)
        r = scaled(row["oec_host_radius_solar"], 100)
        o = scaled(row["oec_reported_a_au"], 1000)
        lines += [
            f"-- {name}",
            f"theorem oecAxis{index} : axisChecks {q} {r} {o} {rounded} := by",
            "  unfold axisChecks",
            "  decide",
            f"#print axioms oecAxis{index}",
        ]
    return lines


def generate(path):
    lines = ["import OpenSpace.Foundations", "",
             "set_option maxRecDepth 4096",
             "set_option maxHeartbeats 2000000",
             "namespace OpenSpace", ""]
    lines += nicer_lines() + contact_lines() + orbit_lines()
    lines += ["", "end OpenSpace", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python3 generate.py OUTPUT.lean")
    generate(Path(sys.argv[1]))
