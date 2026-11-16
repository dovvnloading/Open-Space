"""Run the checks, apply acknowledgements and assemble the ledger."""

from __future__ import annotations

import hashlib
from pathlib import Path
import platform
import sys

from . import catalog, integrity, recompute, replay, trace
from .model import FAIL, INFO, PASS, WARN, Finding, sort_key
from .tree import RepoTree

LEDGER_SCHEMA = "open-space/verification-ledger/1"
DEFAULT_CHECKS = ("catalog", "integrity", "trace", "recompute")
ALL_CHECKS = DEFAULT_CHECKS + ("replay",)


def tool_identity() -> dict:
    """Hash of this package's own source, so a ledger names the exact verifier that produced it."""
    digest = hashlib.sha256()
    folder = Path(__file__).resolve().parent
    files = sorted(folder.glob("*.py"))
    for path in files:
        digest.update(path.name.encode())
        digest.update(b"\0")
        digest.update(path.read_bytes().replace(b"\r\n", b"\n"))
        digest.update(b"\0")
    return {"name": "open_space.verify", "source_sha256": digest.hexdigest(), "files": len(files),
            "python": platform.python_version()}


def _apply_acknowledgements(tree: RepoTree, findings: list[Finding]) -> list[dict]:
    """Downgrade warnings that a claims file acknowledges with a reason. Failures cannot be acknowledged."""
    entries = recompute.acknowledgements(tree)
    used = []
    for entry in entries:
        scope = entry["scope"].rstrip("/") + "/"
        hits = 0
        for finding in findings:
            if finding.severity != WARN or finding.rule != entry["rule"]:
                continue
            if not (finding.path.startswith(scope) or finding.result and f"results/{finding.result}/".startswith(scope)):
                continue
            if entry["key"] and entry["key"] not in finding.key and entry["key"] not in finding.message:
                continue
            finding.severity = INFO
            finding.data = dict(finding.data, acknowledged=entry["reason"], acknowledged_in=entry["source"])
            hits += 1
        used.append(dict(entry, matched=hits))
    return used


def _percent(part: float, whole: float) -> float | None:
    return round(100.0 * part / whole, 1) if whole else None


def scorecards(tree: RepoTree, findings: list[Finding], sections: dict) -> list[dict]:
    """One row per result folder: what was checked and what it showed."""
    labels = catalog.labels_by_result(tree)
    unlisted = {finding.result for finding in findings if finding.rule == "catalog/result-unlisted"}
    traced = sections.get("trace", {}).get("per_result", {})
    replayed = sections.get("replay", {}).get("per_result", {})
    cards: dict[str, dict] = {}
    for name in tree.result_dirs():
        cards[name] = {
            "result": name, "labels": labels.get(name, []), "listed": name not in unlisted,
            "checksums": {"verified": 0, "external": 0, "problems": 0},
            "numbers": None, "claims": {"declared": 0, "reproduced": 0, "failed": 0},
            "replay": replayed.get(name), "fail": 0, "warn": 0,
        }
        numbers = traced.get(name)
        if numbers:
            matched = numbers.get("confirmed", 0) + numbers.get("found", 0)
            cards[name]["numbers"] = {
                "precise": int(numbers.get("precise", 0)), "found": int(matched),
                "confirmed": int(numbers.get("confirmed", 0)), "unmatched": int(numbers.get("unmatched", 0)),
                "mismatch": int(numbers.get("mismatch", 0)),
                "found_percent": _percent(matched, numbers.get("precise", 0)),
                "sensitivity_percent": _percent(matched - numbers.get("expected_accidental_matches", 0.0), matched),
                "evidence_values": int(numbers.get("evidence_values", 0)),
            }
    for finding in findings:
        card = cards.get(finding.result)
        if card is None:
            continue
        if finding.severity == FAIL:
            card["fail"] += 1
        elif finding.severity == WARN:
            card["warn"] += 1
        if finding.rule == "integrity/digests-verified":
            card["checksums"]["verified"] += finding.data["path_bound"] + finding.data["by_content"]
        elif finding.rule == "integrity/external-receipts":
            card["checksums"]["external"] += finding.data["digests"]
        elif finding.check == "integrity" and finding.severity in (FAIL, WARN):
            card["checksums"]["problems"] += 1
        elif finding.check == "recompute":
            card["claims"]["declared"] += 1
            if finding.severity == PASS:
                card["claims"]["reproduced"] += 1
            else:
                card["claims"]["failed"] += 1
    return [cards[name] for name in sorted(cards)]


def run_checks(tree: RepoTree, checks=DEFAULT_CHECKS, only: set[str] | None = None,
               replay_options: dict | None = None, progress=None) -> dict:
    """Run the named checks and return the ledger as a plain dictionary."""
    unknown = [name for name in checks if name not in ALL_CHECKS]
    if unknown:
        raise ValueError(f"unknown checks: {', '.join(unknown)}")
    findings: list[Finding] = []
    sections: dict = {}

    def note(message: str) -> None:
        if progress:
            progress(message)

    if "catalog" in checks:
        note("catalog: links, coverage and entry counts")
        findings.extend(catalog.run(tree))
    if "integrity" in checks:
        note("integrity: recorded checksums against committed bytes")
        produced, extra = integrity.run(tree)
        findings.extend(produced)
        sections["integrity"] = extra["census"]
        sections["external_receipts"] = extra["external"]
    if "trace" in checks:
        note("trace: printed numbers against deposited evidence")
        produced, extra = trace.run(tree, only)
        findings.extend(produced)
        sections["trace"] = extra
    if "recompute" in checks:
        note("recompute: declared claims with independent arithmetic")
        produced, extra = recompute.run(tree, only)
        findings.extend(produced)
        sections["recompute"] = extra
    if "replay" in checks:
        note("replay: documented commands in an offline sandbox")
        produced, extra = replay.run(tree, only, **(replay_options or {}))
        findings.extend(produced)
        sections["replay"] = extra
    if only is not None:
        findings = [finding for finding in findings if not finding.result or finding.result in only]
    acknowledged = _apply_acknowledgements(tree, findings)
    findings.sort(key=sort_key)
    cards = scorecards(tree, findings, sections)
    if only is not None:
        cards = [card for card in cards if card["result"] in only]
    return {
        "schema": LEDGER_SCHEMA,
        "tree": tree.identity(),
        "tool": tool_identity(),
        "checks": [name for name in ALL_CHECKS if name in checks],
        "scope": sorted(only) if only is not None else "all",
        "summary": summarise(findings, sections, cards),
        "acknowledged": acknowledged,
        "results": cards,
        "findings": findings,
        "sections": sections,
    }


def summarise(findings: list[Finding], sections: dict, cards: list[dict]) -> dict:
    counts = {severity: 0 for severity in (FAIL, WARN, INFO, PASS)}
    for finding in findings:
        counts[finding.severity] += 1
    summary: dict = {"findings": counts, "result_folders": len(cards)}
    for finding in findings:
        if finding.rule == "catalog/links-resolve":
            summary["links"] = finding.data
        elif finding.rule in ("catalog/counts-agree", "catalog/count-mismatch"):
            summary["catalog_counts"] = finding.data
        elif finding.rule == "catalog/results-listed":
            summary["catalog_coverage"] = finding.data
    if "integrity" in sections:
        census = sections["integrity"]
        summary["checksums"] = {
            "recorded": census["assertions"], "distinct": census["distinct_digests"],
            "verified_against_committed_bytes": census["path_bound_verified"] + census["content_verified"],
            "mismatched": sum(1 for finding in findings if finding.check == "integrity" and finding.severity == FAIL),
            "external_receipts": census["external"], "external_with_url": census.get("external_with_url", 0),
            "listed_but_not_deposited": census["listed_but_absent"],
        }
    if "trace" in sections:
        totals: dict[str, float] = {}
        for values in sections["trace"]["per_result"].values():
            for key, value in values.items():
                totals[key] = totals.get(key, 0) + value
        matched = totals.get("confirmed", 0) + totals.get("found", 0)
        summary["numbers"] = {
            "precise": int(totals.get("precise", 0)), "found": int(matched),
            "confirmed": int(totals.get("confirmed", 0)), "unmatched": int(totals.get("unmatched", 0)),
            "mismatch": int(totals.get("mismatch", 0)), "last_digit": int(totals.get("last_digit", 0)),
            "found_percent": _percent(matched, totals.get("precise", 0)),
            "sensitivity_percent": _percent(matched - totals.get("expected_accidental_matches", 0.0), matched),
            "identifiers": int(totals.get("identifiers", 0)), "identifiers_found": int(totals.get("identifiers_found", 0)),
            "evidence_values": int(totals.get("evidence_values", 0)),
            "evidence_tables": int(totals.get("evidence_tables", 0)),
            "evidence_documents": int(totals.get("evidence_documents", 0)),
            "catalog": sections["trace"]["catalog"],
        }
    if "recompute" in sections:
        summary["claims"] = sections["recompute"]
    if "replay" in sections:
        summary["replay"] = {key: value for key, value in sections["replay"].items() if key != "per_result"}
    return summary


def to_json(ledger: dict, full: bool = False) -> dict:
    """The ledger as JSON-ready data. Without ``full`` the long per-number match lists are dropped."""
    def slim(finding: Finding) -> dict:
        record = finding.to_json()
        if not full and finding.rule == "trace/numbers-traced" and "data" in record:
            record["data"] = {key: value for key, value in record["data"].items() if key != "matches"}
        return record

    out = {key: value for key, value in ledger.items() if key not in ("findings", "sections")}
    out["findings"] = [slim(finding) for finding in ledger["findings"]
                       if full or finding.severity != PASS or finding.check in ("recompute", "catalog")
                       or finding.rule == "trace/catalog-traced"]
    if full:
        out["external_receipts"] = ledger["sections"].get("external_receipts", {})
    return out


def exit_code(ledger: dict, strict: bool = False) -> int:
    counts = ledger["summary"]["findings"]
    if counts[FAIL] or (strict and counts[WARN]):
        return 1
    return 0


def interpreter() -> str:
    return sys.executable
