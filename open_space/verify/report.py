"""Render a ledger as a Markdown report and as a short console summary."""

from __future__ import annotations

from .model import FAIL, INFO, PASS, WARN, Finding

REPLAY_ORDER = ("reproduced", "ran", "drift", "failed", "needs-network", "needs-dependency", "needs-input",
                "timeout", "script-missing", "template")
REPLAY_WORDS = {
    "reproduced": "reproduced", "ran": "ran, nothing to compare", "drift": "outputs drifted", "failed": "failed",
    "needs-network": "need the network", "needs-dependency": "need a missing package",
    "needs-input": "need inputs that are not deposited", "timeout": "timed out",
    "script-missing": "name a script that does not exist", "template": "contain a placeholder",
}


def _n(value) -> str:
    return f"{value:,}" if isinstance(value, int) else str(value)


def _pct(value) -> str:
    return "n/a" if value is None else f"{value:.1f}%"


def glance(ledger: dict) -> list[tuple[str, str]]:
    """The headline rows: one per check that ran."""
    summary = ledger["summary"]
    rows: list[tuple[str, str]] = []
    if "catalog" in ledger["checks"]:
        parts = []
        counts = summary.get("catalog_counts")
        if counts:
            agree = counts["stated_total"] == counts["entries"] and all(
                counts["stated"].get(label) == found for label, found in counts["labelled_entries"].items()
                if label in counts["stated"])
            parts.append(f"{_n(counts['entries'])} labelled entries; the stated counts "
                         + ("agree with them" if agree else "do NOT agree with them"))
        links = summary.get("links")
        if links:
            parts.append(f"{_n(links['links_checked'])} relative links checked, {_n(links['broken'])} broken")
        coverage = summary.get("catalog_coverage")
        if coverage:
            missing = coverage["folders"] - coverage["listed"]
            parts.append(f"{_n(coverage['listed'])} of {_n(coverage['folders'])} result folders are listed"
                         + (f" ({missing} missing)" if missing else ""))
        rows.append(("Catalog", "; ".join(parts) or "nothing to check"))
    if "checksums" in summary:
        c = summary["checksums"]
        rows.append(("Checksums", f"{_n(c['verified_against_committed_bytes'])} recorded SHA-256 digests match the "
                     f"committed bytes, {_n(c['mismatched'])} do not; {_n(c['external_receipts'])} more describe "
                     f"upstream files that are not in the repository ({_n(c['external_with_url'])} with a source URL)"))
    if "numbers" in summary:
        n = summary["numbers"]
        close = [finding for finding in ledger["findings"] if finding.rule == "trace/mismatch"]
        explained = sum(1 for finding in close if finding.data.get("acknowledged"))
        nearly = f"{_n(len(close) - explained)} sit suspiciously close to a different stored value"
        if explained:
            nearly += f" ({_n(explained)} more do and are explained)"
        rows.append(("Numbers in result READMEs", f"{_n(n['found'])} of {_n(n['precise'])} precise numbers "
                     f"({_pct(n['found_percent'])}) are found in the deposited evidence; {_n(n['unmatched'])} have no "
                     f"counterpart; {nearly}. "
                     f"A mistyped number would have been caught about {_pct(n['sensitivity_percent'])} of the time"))
        cat = n.get("catalog") or {}
        if cat.get("precise"):
            rows.append(("Numbers in the top-level catalog", f"{_n(cat['confirmed'] + cat['found'])} of "
                         f"{_n(cat['precise'])} precise numbers appear in the result they cite"))
    if "claims" in summary:
        c = summary["claims"]
        rows.append(("Claims recomputed independently", f"{_n(c['reproduced'])} of {_n(c['claims'])} declared claims "
                     f"reproduce from committed data, in {_n(c['claim_files'])} result folders"))
    if "replay" in summary:
        r = summary["replay"]
        by = r["by_status"]
        parts = [f"{_n(by[key])} {REPLAY_WORDS[key]}" for key in REPLAY_ORDER if by.get(key)]
        rows.append(("Documented commands replayed offline", f"{_n(r['commands'])} commands: " + ", ".join(parts)))
    return rows


def _location(finding: Finding) -> str:
    return f"`{finding.where()}`" if finding.path else "repository"


def _replay_cell(replay: dict | None) -> str:
    if not replay:
        return ""
    total = sum(replay.values())
    good = replay.get("reproduced", 0)
    bad = replay.get("drift", 0) + replay.get("failed", 0) + replay.get("script-missing", 0)
    cell = f"{good}/{total}"
    return cell + (f" ({bad} to check)" if bad else "")


def markdown(ledger: dict, limit: int = 60) -> str:
    """The full report. Deterministic: the same tree and tool give the same text."""
    tree = ledger["tree"]
    findings: list[Finding] = ledger["findings"]
    lines: list[str] = ["# Verification report", ""]
    commit = tree.get("commit")
    origin = {"commit": f"commit `{commit[:12]}`" if commit else "the committed tree",
              "worktree": "the working tree as it would be committed"
                          + (f" (on top of `{commit[:12]}`)" if commit else ""),
              "directory": "a plain directory", "memory": "an in-memory tree"}[tree["origin"]]
    lines += [
        f"This report was produced by `python -m open_space.verify` from {origin}: "
        f"{_n(tree['files'])} files, {_n(tree['bytes'])} bytes, tree digest `{tree['tree_sha256'][:16]}…`. "
        "Anyone can regenerate it with that command. It does not judge whether a scientific interpretation is "
        "right. Each check answers a narrower question that has exactly one correct answer.", "",
        "## At a glance", "", "| What was checked | Outcome |", "| --- | --- |",
    ]
    lines += [f"| {name} | {outcome} |" for name, outcome in glance(ledger)]
    scope = ledger["scope"]
    if scope != "all":
        lines += ["", f"Scope: only {', '.join(f'`{name}`' for name in scope)}."]
    failures = [finding for finding in findings if finding.severity == FAIL]
    warnings = [finding for finding in findings if finding.severity == WARN]
    lines += ["", "## Needs attention", ""]
    if not failures and not warnings:
        lines += ["Nothing. No check failed and no warning is open.", ""]
    for title, group, lead in (
            ("Failures", failures, "A machine-checkable statement in the repository is false."),
            ("Warnings", warnings, "Not provably wrong, but a person or agent should look.")):
        if not group:
            continue
        lines += [f"### {title} ({len(group)})", "", lead, ""]
        for finding in group[:limit]:
            lines.append(f"- **{finding.rule}** at {_location(finding)}: {finding.message}")
        if len(group) > limit:
            lines.append(f"- … and {len(group) - limit} more in the JSON ledger.")
        lines.append("")
    acknowledged = [finding for finding in findings if finding.data.get("acknowledged")]
    if acknowledged:
        lines += [f"### Acknowledged ({len(acknowledged)})", "",
                  "Warnings that a result's `claims.json` explains. They stay on record with the reason given.", ""]
        for finding in acknowledged[:limit]:
            lines.append(f"- **{finding.rule}** at {_location(finding)}: {finding.message} "
                         f"*Reason: {finding.data['acknowledged']}*")
        lines.append("")

    cards = ledger["results"]
    if cards:
        has_replay = "replay" in ledger["checks"]
        lines += ["## Result by result", "",
                  "Checksums: recorded digests that match committed bytes. Numbers: precise numbers of the folder's "
                  "Markdown found in its deposited evidence. Sensitivity: how often a mistyped number would have been "
                  "caught. Claims: declared claims recomputed independently."
                  + (" Replay: documented commands reproduced offline, out of all documented." if has_replay else ""),
                  ""]
        header = "| Result | Label | Checksums | Numbers found | Sensitivity | Claims |"
        rule = "| --- | --- | ---: | ---: | ---: | ---: |"
        if has_replay:
            header += " Replay |"
            rule += " ---: |"
        lines += [header + " Open |", rule + " ---: |"]
        for card in cards:
            numbers = card["numbers"]
            found = f"{numbers['found']}/{numbers['precise']} ({_pct(numbers['found_percent'])})" \
                if numbers and numbers["precise"] else "none stated"
            sensitivity = _pct(numbers["sensitivity_percent"]) if numbers and numbers["found"] else ""
            claims = card["claims"]
            claim_cell = f"{claims['reproduced']}/{claims['declared']}" if claims["declared"] else ""
            label = ", ".join(card["labels"]) or ("not listed" if not card["listed"] else "")
            open_items = ", ".join(part for part in (
                f"{card['fail']} fail" if card["fail"] else "", f"{card['warn']} warn" if card["warn"] else "") if part)
            row = (f"| [{card['result']}](results/{card['result']}/) | {label} | {card['checksums']['verified']} | "
                   f"{found} | {sensitivity} | {claim_cell} |")
            if has_replay:
                row += f" {_replay_cell(card['replay'])} |"
            lines.append(row + f" {open_items} |")
        lines.append("")

    lines += [
        "## What each check shows, and what it does not", "",
        "- **Catalog.** Every relative link in every Markdown file resolves, every result folder is reachable from "
        "the top-level README, and the entry counts stated there equal the labelled entries. It says nothing about "
        "whether an entry's wording is fair.",
        "- **Checksums.** Each SHA-256 recorded for a committed file equals the hash of the committed bytes, read "
        "from the Git object store so that line-ending conversion cannot interfere. Digests of upstream downloads "
        "are counted; they can only be checked by fetching those files again.",
        "- **Numbers.** Each number a README prints, with the precision it prints, is searched for in the folder's "
        "JSON, tables and scripts, including facts computed from tables such as row counts and extremes. A match "
        "means the prose and the deposited evidence agree on that number. It does not mean the number is the right "
        "measurement of the sky. An unmatched number is not an error: it has no machine-readable counterpart yet.",
        "- **Claims.** Where a folder declares claims in `claims.json`, each one is recomputed from committed data "
        "with this package's own standard-library arithmetic, never the code that produced it, and must round to "
        "exactly what the README prints.",
        "- **Replay.** Each `python` command a README documents is run in a sandbox holding the committed tree, with "
        "the network disabled, and the files it rewrites are compared with the committed ones.", "",
        "## Reproduce this report", "",
        "```sh",
        "python -m open_space.verify --report VERIFICATION.md --json verification-ledger.json",
        "```", "",
        f"Verifier source digest `{ledger['tool']['source_sha256'][:16]}…`. The same tree and the same verifier "
        "give the same report on any machine.", "",
    ]
    return "\n".join(lines)


def console(ledger: dict, limit: int = 25) -> str:
    """A short plain-text summary for the terminal."""
    tree = ledger["tree"]
    findings: list[Finding] = ledger["findings"]
    counts = ledger["summary"]["findings"]
    where = tree["commit"][:12] if tree.get("commit") else tree["origin"]
    out = [f"open_space.verify  {tree['origin']} {where}  {_n(tree['files'])} files  tree {tree['tree_sha256'][:12]}", ""]
    width = max((len(name) for name, _ in glance(ledger)), default=0)
    for name, outcome in glance(ledger):
        out.append(f"  {name.ljust(width)}  {outcome}")
    out.append("")
    for severity, title in ((FAIL, "FAIL"), (WARN, "WARN")):
        group = [finding for finding in findings if finding.severity == severity]
        for finding in group[:limit]:
            out.append(f"  {title}  {finding.rule}  {finding.where()}")
            out.append(f"        {finding.message}")
        if len(group) > limit:
            out.append(f"  {title}  … and {len(group) - limit} more")
    out.append("")
    out.append(f"  {counts[FAIL]} failed, {counts[WARN]} warnings, {counts[INFO]} notes, {counts[PASS]} checks passed")
    return "\n".join(out)
