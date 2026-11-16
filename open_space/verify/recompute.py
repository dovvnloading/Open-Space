"""Recompute declared claims from committed data with independent arithmetic.

A result folder may carry a ``claims.json``. Each entry ties one number printed
in the README to an expression that recomputes it from the folder's own tables
and JSON, using only the estimators in :mod:`open_space.verify.stats`::

    {
      "schema": "open-space/claims/1",
      "claims": [
        {
          "id": "hard-band-correlation",
          "printed": "0.96857",
          "let": {"t": "table('profiles.csv')"},
          "compute": "pearson(t['stix_50_84'], t['gbm_channels_2_4_rate'])"
        }
      ]
    }

The claim holds when the recomputed value rounds to exactly what is printed.
The printed text must still be present in the document, so a README edit that
changes a number cannot leave a stale claim passing.
"""

from __future__ import annotations

import json
import math
import posixpath

from . import numbers
from .expr import Evaluator, ExprError, Vec
from .model import FAIL, PASS, Finding
from .tree import RepoTree

CHECK = "recompute"
SCHEMA = "open-space/claims/1"
FILENAME = "claims.json"


def claim_files(tree: RepoTree) -> list[str]:
    return [path for path in tree.paths if posixpath.basename(path) == FILENAME]


def load(tree: RepoTree, path: str) -> tuple[dict | None, str | None]:
    """Parse a claims file. Returns (document, error)."""
    try:
        document = json.loads(tree.text(path))
    except ValueError as error:
        return None, f"not valid JSON: {error}"
    if not isinstance(document, dict):
        return None, "top level must be an object"
    if document.get("schema") != SCHEMA:
        return None, f"schema must be {SCHEMA!r}"
    for field, kind in (("claims", list), ("acknowledged", list), ("let", dict)):
        if field in document and not isinstance(document[field], kind):
            return None, f"{field!r} must be a {kind.__name__}"
    return document, None


def acknowledgements(tree: RepoTree) -> list[dict]:
    """Every acknowledged finding in the tree: {'scope', 'rule', 'key', 'reason', 'source'}."""
    found = []
    for path in claim_files(tree):
        document, error = load(tree, path)
        if error or not document:
            continue
        scope = posixpath.dirname(path)
        for entry in document.get("acknowledged", []):
            if isinstance(entry, dict) and isinstance(entry.get("rule"), str) and isinstance(entry.get("reason"), str) \
                    and entry["reason"].strip():
                found.append({"scope": scope, "rule": entry["rule"], "key": str(entry.get("key", "")),
                              "reason": entry["reason"].strip(), "source": path})
    return found


def _printed_number(text: str) -> numbers.Number | None:
    found = [number for number in numbers.extract(text) if number.kind in ("value", "year", "identifier", "date")]
    return found[0] if len(found) == 1 else None


def _is_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _as_scalar(value):
    if isinstance(value, Vec):
        if len(value) != 1:
            raise ExprError(f"expression produced a column of {len(value)} values, not one number")
        value = value.data[0]
    if isinstance(value, (list, tuple)):
        if len(value) != 1:
            raise ExprError(f"expression produced {len(value)} values, not one")
        value = value[0]
    return value


def _agrees(computed: float, printed: numbers.Number, tolerance) -> tuple[bool, str]:
    """Compare with the printed value. Default: the computed value must round to what is printed."""
    difference = abs(computed - printed.value)
    if tolerance is None:
        candidates = [(printed.value, printed.half)] + [(value, ulp / 2) for value, ulp, _ in printed.readings]
        ok = any(numbers.round_matches(computed, value, half) for value, half in candidates)
        return ok, "rounds to the printed value"
    if not isinstance(tolerance, dict) or len(tolerance) != 1:
        raise ExprError("tolerance must be one of {'abs': x}, {'rel': x} or {'ulps': n}")
    (kind, amount), = tolerance.items()
    if not isinstance(amount, (int, float)) or isinstance(amount, bool) or amount < 0:
        raise ExprError("tolerance amount must be a non-negative number")
    if kind == "abs":
        return difference <= amount, f"within {amount:g}"
    if kind == "rel":
        return difference <= amount * abs(printed.value), f"within a relative {amount:g}"
    if kind == "ulps":
        return difference <= amount * printed.ulp * (1 + 1e-9), f"within {amount:g} units in the last printed place"
    raise ExprError(f"unknown tolerance kind {kind!r}")


def check_file(tree: RepoTree, path: str) -> list[Finding]:
    """Evaluate every claim of one claims file."""
    document, error = load(tree, path)
    if error:
        return [Finding(CHECK, "recompute/spec-invalid", FAIL, path, f"claims file is unusable: {error}", key="file")]
    folder = posixpath.dirname(path)
    findings: list[Finding] = []
    evaluator = Evaluator(tree.read, folder)
    shared_error = None
    for name, source in (document.get("let") or {}).items():
        try:
            evaluator.define(name, source)
        except ExprError as failure:
            shared_error = f"shared definition {name!r} failed: {failure}"
            break
    seen: set[str] = set()
    for position, claim in enumerate(document.get("claims", [])):
        identifier = claim.get("id") if isinstance(claim, dict) else None
        if not isinstance(identifier, str) or not identifier or identifier in seen:
            findings.append(Finding(CHECK, "recompute/spec-invalid", FAIL, path,
                                    f"claim #{position + 1} needs a unique string 'id'", key=f"#{position}"))
            continue
        seen.add(identifier)
        findings.append(_check_claim(tree, path, folder, evaluator, claim, identifier, shared_error))
    return findings


def _check_claim(tree: RepoTree, path: str, folder: str, shared: Evaluator, claim: dict, identifier: str,
                 shared_error: str | None) -> Finding:
    def fail(rule: str, message: str, **data) -> Finding:
        return Finding(CHECK, rule, FAIL, path, f"{identifier}: {message}", key=identifier,
                       data={"claim": identifier, **data})

    source = claim.get("compute")
    if not isinstance(source, str) or not source.strip():
        return fail("recompute/spec-invalid", "needs a 'compute' expression")
    has_printed = "printed" in claim
    has_expect = "expect" in claim
    if has_printed == has_expect:
        return fail("recompute/spec-invalid", "needs exactly one of 'printed' (text in the document) or 'expect'")
    document_path = posixpath.normpath(posixpath.join(folder, claim.get("document", "README.md")))
    printed = None
    line = None
    if has_printed:
        text = claim["printed"]
        if not isinstance(text, str) or not text.strip():
            return fail("recompute/spec-invalid", "'printed' must be the text as it appears in the document")
        if document_path not in tree:
            return fail("recompute/claim-not-in-document", f"document {document_path} does not exist")
        body = tree.text(document_path)
        at = body.find(text)
        if at < 0:
            return fail("recompute/claim-not-in-document",
                        f"the text {text!r} is no longer in {document_path}; the claim is stale", printed=text)
        line = body.count("\n", 0, at) + 1
        printed = _printed_number(text)
        if printed is None:
            return fail("recompute/spec-invalid", f"'printed' text {text!r} must contain exactly one number")
    if shared_error:
        return fail("recompute/error", shared_error)
    evaluator = Evaluator(tree.read, folder)
    evaluator.names.update(shared.names)
    evaluator._tables, evaluator._documents = shared._tables, shared._documents      # reuse parsed files
    try:
        for name, definition in (claim.get("let") or {}).items():
            evaluator.define(name, definition)
        value = _as_scalar(evaluator.evaluate(source))
    except ExprError as failure:
        return fail("recompute/error", f"could not be recomputed: {failure}", compute=source)
    statement = claim.get("statement", "")
    if has_expect:
        expected = claim["expect"]
        if isinstance(value, Vec):
            value = value.data
        numeric = _is_number(expected) and _is_number(value)
        if numeric and (isinstance(expected, float) or isinstance(value, float) or "tolerance" in claim):
            tolerance = claim.get("tolerance") or {"rel": 1e-12}
            if not isinstance(tolerance, dict) or len(tolerance) != 1 or next(iter(tolerance)) not in ("abs", "rel"):
                return fail("recompute/spec-invalid", "'expect' takes a tolerance of {'abs': x} or {'rel': x}")
            (kind, amount), = tolerance.items()
            if not _is_number(amount) or amount < 0:
                return fail("recompute/spec-invalid", "tolerance amount must be a non-negative number")
            agrees = math.isfinite(value) and abs(value - expected) <= amount * (abs(expected) if kind == "rel" else 1.0)
            how = f"within {kind} {amount:g} of the expected value"
        else:
            # Whole numbers, truth values, text and lists must match exactly, and True is not 1.
            agrees = value == expected and isinstance(value, bool) == isinstance(expected, bool)
            how = "equals the expected value"
        if not agrees:
            return fail("recompute/disagrees", f"recomputed {value!r}, expected {expected!r}", computed=repr(value),
                        expected=expected, compute=source)
        return Finding(CHECK, "recompute/reproduced", PASS, path, f"{identifier}: recomputed {value!r}, {how}",
                       key=identifier, data={"claim": identifier, "computed": value, "expected": expected,
                                             "statement": statement, "compute": source})
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return fail("recompute/error", f"expression produced {value!r}, not a finite number", compute=source)
    try:
        agrees, how = _agrees(float(value), printed, claim.get("tolerance"))
    except ExprError as failure:
        return fail("recompute/spec-invalid", str(failure))
    data = {"claim": identifier, "printed": claim["printed"], "computed": float(value), "document": document_path,
            "line": line, "statement": statement, "compute": source}
    if not agrees:
        return Finding(CHECK, "recompute/disagrees", FAIL, path,
                       f"{identifier}: the document prints {claim['printed']}, independent recomputation gives "
                       f"{float(value):.12g}", key=identifier, data=data)
    return Finding(CHECK, "recompute/reproduced", PASS, path,
                   f"{identifier}: printed {claim['printed']}, recomputed {float(value):.12g} ({how})",
                   key=identifier, data=data)


def run(tree: RepoTree, only: set[str] | None = None) -> tuple[list[Finding], dict]:
    findings: list[Finding] = []
    for path in claim_files(tree):
        parts = path.split("/")
        if only is not None and not (len(parts) >= 3 and parts[0] == "results" and parts[1] in only):
            continue
        findings.extend(check_file(tree, path))
    summary = {
        "claim_files": len({finding.path for finding in findings}),
        "claims": sum(1 for finding in findings if finding.rule != "recompute/spec-invalid" or finding.key != "file"),
        "reproduced": sum(1 for finding in findings if finding.severity == PASS),
        "failed": sum(1 for finding in findings if finding.severity == FAIL),
    }
    return findings, summary
