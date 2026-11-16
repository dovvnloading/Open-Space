"""Findings: the one record type every check emits."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib

FAIL = "fail"   # a machine-checkable assertion in the repository is false
WARN = "warn"   # ambiguous or suspicious; a person or agent should look
INFO = "info"   # a gap in what can be checked mechanically, not an error
PASS = "pass"   # an assertion that was checked and holds

SEVERITIES = (FAIL, WARN, INFO, PASS)


def result_of(path: str) -> str:
    """Return the result folder a path belongs to, or '' for repository-level files."""
    parts = path.split("/")
    if len(parts) >= 3 and parts[0] == "results":
        return parts[1]
    return ""


@dataclass
class Finding:
    """One verdict about one assertion, with enough data to re-check it by hand."""

    check: str                      # family: catalog, integrity, trace, recompute, replay
    rule: str                       # stable rule name, for example "integrity/digest-mismatch"
    severity: str
    path: str                       # file that carries the assertion
    message: str
    line: int | None = None
    key: str = ""                   # distinguishes findings of one rule in one file without using line numbers
    data: dict = field(default_factory=dict)
    result: str = ""

    def __post_init__(self) -> None:
        if self.severity not in SEVERITIES:
            raise ValueError(f"unknown severity {self.severity!r}")
        if not self.result:
            self.result = result_of(self.path)

    @property
    def id(self) -> str:
        """Identity that survives unrelated edits: rule, file and key, never the line number."""
        raw = "\0".join((self.rule, self.path, self.key)).encode("utf-8", "surrogateescape")
        return hashlib.sha1(raw).hexdigest()[:12]

    def to_json(self) -> dict:
        out = {
            "id": self.id,
            "check": self.check,
            "rule": self.rule,
            "severity": self.severity,
            "result": self.result,
            "path": self.path,
            "message": self.message,
        }
        if self.line is not None:
            out["line"] = self.line
        if self.data:
            out["data"] = self.data
        return out

    def where(self) -> str:
        return self.path if self.line is None else f"{self.path}:{self.line}"


def sort_key(finding: Finding) -> tuple:
    return (SEVERITIES.index(finding.severity), finding.check, finding.rule, finding.path,
            finding.line or 0, finding.key)
