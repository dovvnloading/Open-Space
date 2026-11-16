"""Trace every number a README prints to a value in the deposited evidence.

The question for each printed number is narrow: does the folder hold a
machine-readable value that rounds to exactly what the prose says?

A match is only evidence when a *wrong* number would not also have matched. The
check therefore measures, for every match, how crowded the evidence is around
the printed value: the fraction of neighbouring rounding cells that hold some
value. That fraction is the probability that a number of the same precision,
picked nearby at random, would have matched by accident. It is added up over
every place that was searched before the match.

Outcomes for precise numbers (three or more significant digits):

* confirmed: found, and a wrong number would have matched at most 5% of the time;
* found: a value rounds to it, but the evidence is crowded there, so the match
  alone proves less;
* unmatched: nothing in the folder rounds to it.

One minus the average chance over the matched numbers is the sensitivity of the
check: how often a mistyped number would have been caught.

A mismatch alarm is raised only when it could hardly be a coincidence: the
printed number sits a few last-place units from an evidence value, or differs
from one by a single digit or a swapped pair, and the evidence is so sparse at
that precision that an unrelated number would land there less than once in a
thousand times.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass
import math

from . import numbers as num
from .catalog import ROOT_README, catalog_entries
from .evidence import Evidence, Pool
from .model import INFO, PASS, WARN, Finding
from .numbers import Number
from .tree import RepoTree

CHECK = "trace"
ALPHA = 0.05            # largest accepted chance that a wrong number would have matched
CELLS = 50              # rounding cells inspected on each side of a printed value
NEAR_ULPS = 10          # how many last-place units away a transcription slip may sit
NEAR_ALPHA = 0.001      # largest accepted chance that a mismatch alarm is a coincidence
_SQRT10 = math.sqrt(10.0)
_LN10 = math.log(10.0)

# Unit conversions tried for precise numbers after every literal reading has failed.
SCALES = (
    (100.0, "×100"), (0.01, "÷100"), (1e3, "×1000"), (1e-3, "÷1000"), (1e6, "×10⁶"), (1e-6, "÷10⁶"),
    (60.0, "×60"), (1 / 60.0, "÷60"), (3600.0, "×3600"), (1 / 3600.0, "÷3600"),
    (24.0, "×24"), (1 / 24.0, "÷24"), (1440.0, "×1440"), (1 / 1440.0, "÷1440"),
    (86400.0, "×86400"), (1 / 86400.0, "÷86400"), (365.25, "×365.25"), (1 / 365.25, "÷365.25"),
)

EXACT, LOOSE, CONVERTED = 0, 1, 2
LOOSE_NOTE = "last printed digit differs by rounding or truncation"


@dataclass
class Trace:
    """What the search found for one printed number."""

    number: Number
    status: str                  # confirmed, found, mismatch or unmatched
    chance: float = 0.0          # probability a wrong number would have matched by this point
    pool: str = ""
    where: str = ""
    found: float | None = None   # the evidence value
    how: str = ""                # reading that matched: "as printed", "sign dropped in prose", "evidence ×100" ...
    note: str = ""               # LOOSE_NOTE, or what kind of mismatch


def _slack(target: float, half: float) -> float:
    return half * 1e-9 + abs(target) * 4 * math.ulp(1.0)


def _candidates(values: list[float], target: float, half: float) -> range:
    """Indices of the values inside ``target ± half``."""
    slack = _slack(target, half)
    low = bisect.bisect_left(values, target - half - slack)
    high = bisect.bisect_right(values, target + half + slack)
    return range(low, high)


def _closest(values: list[float], target: float, indices) -> int | None:
    best = None
    for index in indices:
        if best is None or abs(values[index] - target) < abs(values[best] - target):
            best = index
    return best


def occupancy(values: list[float], target: float, half: float, cells: int = CELLS) -> float:
    """Fraction of rounding cells around ``target`` that contain a value.

    Cells are the intervals a number of this precision could round into, ``cells`` of
    them on each side, not counting the target's own.
    """
    if not values or half <= 0:
        return 0.0
    width = 2 * half
    slack = _slack(target, half)
    lowest = target - (cells + 0.5) * width
    highest = target + (cells + 0.5) * width
    start = bisect.bisect_left(values, lowest)
    stop = bisect.bisect_right(values, highest)
    if stop - start > 4 * cells:
        occupied = 0                                   # dense: probe each cell with one bisection
        for step in range(1, cells + 1):
            for side in (-1.0, 1.0):
                low = target + side * step * width - half
                at = bisect.bisect_left(values, low + slack, start, stop)
                if at < stop and values[at] < low + width - slack:
                    occupied += 1
        return occupied / (2 * cells)
    seen = set()                                       # sparse: bin the few values that are there
    for index in range(start, stop):
        cell = math.floor((values[index] - target) / width + 0.5)
        if cell != 0 and -cells <= cell <= cells:
            seen.add(cell)
    return len(seen) / (2 * cells)


def decade_chance(values: list[float], target: float, ulp: float, cells: float) -> float:
    """Chance that an unrelated number lands within ``cells`` last-place units of some value.

    Uses the number of evidence values within half a decade either side of the target,
    spread evenly in logarithm. Unlike ``occupancy`` this resolves very small chances.
    """
    if target == 0 or not values:
        return 1.0
    low, high = sorted((target / _SQRT10, target * _SQRT10))
    count = bisect.bisect_right(values, high) - bisect.bisect_left(values, low)
    return min(1.0, cells * max(count, 1) * ulp / (abs(target) * _LN10))


def _places(number: Number) -> int:
    return -int(round(math.log10(number.ulp))) if number.ulp > 0 else 0


def _readings(number: Number) -> list[tuple[float, float, str]]:
    """Every value the printed text may stand for, most literal first."""
    readings = [(number.value, number.half, "as printed")]
    readings += [(value, ulp / 2, why) for value, ulp, why in number.readings]
    if not number.signed:
        readings += [(-value, half, "sign dropped in prose" if why == "as printed" else why + ", sign dropped in prose")
                     for value, half, why in list(readings) if value > 0]
    return readings


def trace_number(number: Number, evidence: Evidence, extra_pools: tuple[Pool, ...] = (),
                 describe_extra=None) -> Trace:
    """Search the pools in order and return the first match with its accumulated chance."""
    pools = [pool for pool in tuple(extra_pools) + evidence.pools if len(pool)]
    plan: list[tuple[Pool, float, float, str, int]] = []         # pool, target, half-width, reading, stage
    readings = _readings(number)
    for pool in pools:
        for value, half, why in readings:
            plan.append((pool, value, half, why, EXACT))
    if number.precise:
        convertible = [pool for pool in pools if pool is not evidence.code]
        literal = [(number.value, "as printed")]
        if not number.signed and number.value > 0:
            literal.append((-number.value, "sign dropped in prose"))
        for pool in convertible:
            for value, why in literal:
                plan.append((pool, value, number.ulp, why, LOOSE))
        for scale, label in SCALES:
            for pool in convertible:
                for value, why in literal:
                    suffix = "" if why == "as printed" else ", sign dropped in prose"
                    plan.append((pool, value / scale, number.half / abs(scale), f"evidence {label}{suffix}", CONVERTED))
    claim_places = _places(number)
    for position, (pool, target, half, why, stage) in enumerate(plan):
        indices = _candidates(pool.values, target, half)
        if not indices:
            continue
        if stage == LOOSE:
            # Only a value written with more digits than the prose can be a rounding variant of it.
            indices = [index for index in indices
                       if pool.places(index) is None or pool.places(index) > claim_places]
            if not indices:
                continue
        index = _closest(pool.values, target, indices)
        chance = 0.0
        ceiling = ALPHA if stage == CONVERTED else 1.0
        for earlier_pool, earlier_target, earlier_half, _, _ in plan[:position + 1]:
            chance += occupancy(earlier_pool.values, earlier_target, earlier_half)
            if chance > ceiling:
                break
        chance = min(1.0, chance)
        if stage == CONVERTED and chance > ALPHA:
            # A unit conversion is a guess. It only counts when nothing else could explain the match.
            return Trace(number, "unmatched")
        if pool in evidence.pools:
            where = evidence.describe(pool, index)
        else:
            where = describe_extra(pool, index) if describe_extra else pool.name
        return Trace(number, "confirmed" if chance <= ALPHA else "found", chance, pool.name, where,
                     pool.values[index], why, LOOSE_NOTE if stage == LOOSE else "")
    return Trace(number, "unmatched")


def _mantissa(text: str) -> float:
    cleaned = text.replace(",", "").replace(num.MINUS, "-").replace("+", "").strip()
    try:
        return float(cleaned)
    except ValueError:
        return 0.0


def _digit_variants(number: Number) -> list[tuple[float, str]]:
    """Values one typing slip away: a single digit changed, or two neighbours swapped."""
    text = number.text
    mantissa_end = len(text)
    for marker in ("×", "\\times", "\\cdot", "·", "e", "E", "x", "X", "*"):
        at = text.find(marker)
        if at > 0:
            mantissa_end = min(mantissa_end, at)
    head = text[:mantissa_end]
    positions = [index for index, char in enumerate(head) if char.isdigit()]
    base = _mantissa(head)
    if not base or len(positions) < 3:
        return []
    scale = number.value / base
    variants: list[tuple[float, str]] = []
    for place, index in enumerate(positions):
        for digit in "0123456789":
            if digit == head[index] or (place == 0 and digit == "0"):
                continue
            variants.append((_mantissa(head[:index] + digit + head[index + 1:]) * scale, "one digit differs"))
    for first, second in zip(positions, positions[1:]):
        if second == first + 1 and head[first] != head[second]:
            swapped = head[:first] + head[second] + head[first] + head[second + 1:]
            variants.append((_mantissa(swapped) * scale, "two neighbouring digits are swapped"))
    unique: list[tuple[float, str]] = []
    seen = {round(number.value / number.ulp)}
    for value, why in variants:
        key = round(value / number.ulp)
        if key not in seen and math.isfinite(value):
            seen.add(key)
            unique.append((value, why))
    return unique


def mismatch(number: Number, evidence: Evidence) -> Trace | None:
    """For an unmatched number: is there a stored value it almost certainly was meant to be?

    Only values stored in JSON or small tables are considered, and only when the evidence
    is sparse enough at this precision that the resemblance is not a coincidence.
    """
    pool = evidence.direct
    if not number.precise or not len(pool) or number.value == 0:
        return None
    values = pool.values
    claim_places = _places(number)
    targets = [(number.value, "as printed")]
    if not number.signed and number.value > 0:
        targets.append((-number.value, "sign dropped in prose"))

    for target, how in targets:
        index = pool.nearest(target)
        if index is None:
            continue
        gap = abs(values[index] - target) / number.ulp
        written = pool.places(index)
        if 1.0 - 1e-9 <= gap <= NEAR_ULPS and (written is None or written >= claim_places):
            chance = decade_chance(values, target, number.ulp, 2 * NEAR_ULPS)
            if chance <= NEAR_ALPHA:
                return Trace(number, "mismatch", chance, pool.name, evidence.describe(pool, index), values[index],
                             how, f"differs by {gap:.1f} units in the last printed place")
    variants = _digit_variants(number)
    hits = []
    for value, why in variants:
        index = _closest(values, value, _candidates(values, value, number.half))
        if index is not None and pool.places(index) == claim_places:
            hits.append((value, why, index))
    if len(hits) == 1:
        chance = min(1.0, sum(decade_chance(values, value, number.ulp, 1.0) for value, _ in variants))
        if chance <= NEAR_ALPHA:
            value, why, index = hits[0]
            return Trace(number, "mismatch", chance, pool.name, evidence.describe(pool, index), values[index],
                         "as printed", why)
    return None


def _prose_pool(tree: RepoTree, folder: str) -> tuple[Pool, list[str]]:
    """Numbers printed in a result's own Markdown, used to check the top-level catalog against it."""
    pool = Pool("result prose")
    sources: list[str] = []
    for path in tree.under(folder):
        if not path.lower().endswith(".md"):
            continue
        sources.append(path)
        for number in num.extract(tree.text(path)):
            if number.kind != "value":
                continue
            pool.add(number.value, len(sources) - 1, f"line {number.line}", _places(number))
            for value, ulp, _ in number.readings:
                pool.add(value, len(sources) - 1, f"line {number.line}", -int(round(math.log10(ulp))) if ulp > 0 else 0)
    pool.freeze()
    return pool, sources


def trace_document(text: str, evidence: Evidence, extra_pools: tuple[Pool, ...] = (), describe_extra=None,
                   ) -> tuple[list[Number], list[Trace], list[tuple[Number, str | None]]]:
    """Trace every number of one Markdown document against one evidence set."""
    numbers = num.extract(text)
    traces: list[Trace] = []
    identifiers: list[tuple[Number, str | None]] = []
    for number in numbers:
        if number.kind == "identifier":
            digits = "".join(char for char in number.text if char.isdigit())
            if len(digits) >= 5:
                identifiers.append((number, evidence.find_literal(digits)))
            continue
        if number.kind != "value":
            continue
        trace = trace_number(number, evidence, extra_pools, describe_extra)
        if trace.status == "unmatched" and number.precise:
            trace = mismatch(number, evidence) or trace
        traces.append(trace)
    return numbers, traces, identifiers


def summarise(traces: list[Trace]) -> dict:
    """Counts for one document. Sensitivity follows from them: 1 - accidental / matched."""
    precise = [trace for trace in traces if trace.number.precise]
    weak = [trace for trace in traces if trace.number.kind == "value" and not trace.number.precise]
    matched = [trace for trace in precise if trace.status in ("confirmed", "found")]
    return {
        "precise": len(precise),
        "confirmed": sum(1 for trace in precise if trace.status == "confirmed"),
        "found": sum(1 for trace in precise if trace.status == "found"),
        "mismatch": sum(1 for trace in precise if trace.status == "mismatch"),
        "unmatched": sum(1 for trace in precise if trace.status == "unmatched"),
        "last_digit": sum(1 for trace in matched if trace.note == LOOSE_NOTE),
        # Sum of the chances: how many of the matches are expected to be accidents.
        "expected_accidental_matches": round(sum(trace.chance for trace in matched), 3),
        "low_precision": len(weak),
        "low_precision_found": sum(1 for trace in weak if trace.status in ("confirmed", "found")),
    }


def _document_findings(path: str, numbers: list[Number], traces: list[Trace],
                       identifiers: list[tuple[Number, str | None]], result: str) -> list[Finding]:
    findings: list[Finding] = []
    summary = summarise(traces)
    summary["identifiers"] = len(identifiers)
    summary["identifiers_found"] = sum(1 for _, where in identifiers if where)
    summary["not_checked"] = {kind: sum(1 for number in numbers if number.kind == kind)
                              for kind in ("year", "date", "reference")}
    occurrence: dict[str, int] = {}
    for trace in traces:
        if trace.status != "mismatch":
            continue
        number = trace.number
        occurrence[number.text] = occurrence.get(number.text, 0) + 1
        findings.append(Finding(
            CHECK, "trace/mismatch", WARN, path,
            f"prints {number.text}; the evidence holds {trace.found!r} at {trace.where} ({trace.note})",
            line=number.line, key=f"{number.text}#{occurrence[number.text]}", result=result,
            data={"printed": number.text, "evidence_value": trace.found, "evidence": trace.where,
                  "difference": trace.note, "chance_of_coincidence": float(f"{trace.chance:.2g}"),
                  "context": number.context}))
    loose = [trace for trace in traces if trace.note == LOOSE_NOTE and trace.number.precise]
    if loose:
        findings.append(Finding(
            CHECK, "trace/last-digit", INFO, path,
            f"{len(loose)} numbers match the evidence only after allowing one unit in the last printed place",
            key="last-digit", result=result,
            data={"numbers": [{"printed": trace.number.text, "line": trace.number.line, "evidence_value": trace.found,
                               "evidence": trace.where} for trace in loose[:200]], "total": len(loose)}))
    unmatched = [trace.number for trace in traces if trace.status == "unmatched" and trace.number.precise]
    if unmatched:
        findings.append(Finding(
            CHECK, "trace/unmatched", INFO, path,
            f"{len(unmatched)} precise numbers have no counterpart in the deposited evidence", key="unmatched",
            result=result, data={"numbers": [{"printed": number.text, "line": number.line, "context": number.context}
                                             for number in unmatched[:400]], "total": len(unmatched)}))
    missing = [number for number, where in identifiers if not where]
    if missing:
        findings.append(Finding(
            CHECK, "trace/identifier-unmatched", INFO, path,
            f"{len(missing)} identifiers do not appear in any evidence file", key="identifiers", result=result,
            data={"identifiers": [{"printed": number.text, "line": number.line, "context": number.context}
                                  for number in missing[:200]], "total": len(missing)}))
    matched = [trace for trace in traces if trace.status in ("confirmed", "found") and trace.number.precise]
    summary["matches"] = [
        {"printed": trace.number.text, "line": trace.number.line, "evidence": trace.where, "how": trace.how,
         "status": trace.status, "chance": round(trace.chance, 4)} for trace in matched[:600]]
    findings.append(Finding(
        CHECK, "trace/numbers-traced", PASS, path,
        f"{summary['confirmed'] + summary['found']} of {summary['precise']} precise numbers found in the evidence "
        f"({summary['confirmed']} unambiguously)", key="traced", result=result, data=summary))
    return findings


def run(tree: RepoTree, only: set[str] | None = None) -> tuple[list[Finding], dict]:
    """Trace every result folder's Markdown, then the top-level catalog against the folders it cites."""
    findings: list[Finding] = []
    cache: dict[str, Evidence] = {}
    prose_cache: dict[str, tuple[Pool, list[str]]] = {}
    per_result: dict[str, dict] = {}

    def evidence_for(name: str) -> Evidence:
        if name not in cache:
            cache[name] = Evidence.build(tree, f"results/{name}")
        return cache[name]

    def prose_for(name: str) -> tuple[Pool, list[str]]:
        if name not in prose_cache:
            prose_cache[name] = _prose_pool(tree, f"results/{name}")
        return prose_cache[name]

    for folder in tree.result_dirs():
        if only is not None and folder not in only:
            continue
        documents = [path for path in tree.under(f"results/{folder}") if path.lower().endswith(".md")]
        if not documents:
            continue
        evidence = evidence_for(folder)
        totals: dict[str, float] = {}
        for path in documents:
            numbers, traces, identifiers = trace_document(tree.text(path), evidence)
            produced = _document_findings(path, numbers, traces, identifiers, folder)
            findings.extend(produced)
            for key, value in produced[-1].data.items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    totals[key] = totals.get(key, 0) + value
        totals["expected_accidental_matches"] = round(totals.get("expected_accidental_matches", 0.0), 3)
        totals["evidence_values"] = sum(len(pool) for pool in evidence.pools)
        totals["evidence_tables"] = evidence.tables
        totals["evidence_documents"] = evidence.documents
        totals["evidence_scripts"] = evidence.scripts
        per_result[folder] = totals

    catalog = {"entries": 0, "precise": 0, "confirmed": 0, "found": 0, "unmatched": 0}
    if ROOT_README in tree:
        for entry in catalog_entries(tree):
            folders = [name for name in entry["results"] if only is None or name in only]
            if not folders:
                continue
            catalog["entries"] += 1
            unmatched: list[Number] = []
            for number in num.extract(entry["text"]):
                if not number.precise:
                    continue
                catalog["precise"] += 1
                best: Trace | None = None
                for name in folders:
                    prose, sources = prose_for(name)
                    trace = trace_number(
                        number, evidence_for(name), (prose,),
                        lambda pool, index, sources=sources: f"{sources[pool.source(index)]} {pool.detail(index)}")
                    if trace.status in ("confirmed", "found"):
                        best = trace
                        break
                if best is None:
                    catalog["unmatched"] += 1
                    unmatched.append(number)
                else:
                    catalog[best.status] += 1
            if unmatched:
                printed = ", ".join(number.text for number in unmatched)
                findings.append(Finding(
                    CHECK, "trace/catalog-number-unmatched", WARN, ROOT_README,
                    f"catalog entry prints {printed}, which appears neither in the prose nor in the evidence "
                    f"of results/{folders[0]}", line=entry["line"], key=f"{folders[0]}:{printed}", result=folders[0],
                    data={"numbers": [number.text for number in unmatched], "entry": entry["text"][:300]}))
        findings.append(Finding(
            CHECK, "trace/catalog-traced", PASS, ROOT_README,
            f"{catalog['confirmed'] + catalog['found']} of {catalog['precise']} precise numbers in catalog entries "
            "appear in the result they cite", key="catalog", data=catalog))
    return findings, {"per_result": per_result, "catalog": catalog}
