"""Catalog checks: does the public README agree with the files behind it?

Three questions, all answerable without judgement:

* does every relative link in every Markdown file lead to a file that exists;
* is every result folder reachable from the top-level catalog;
* do the entry counts stated in the README equal the labelled entries below it.
"""

from __future__ import annotations

import posixpath
import re
from urllib.parse import unquote

from . import markdown
from .model import FAIL, INFO, PASS, WARN, Finding
from .tree import RepoTree

CHECK = "catalog"
ROOT_README = "README.md"

# Label printed on each catalog bullet, and the words the summary sentence uses for it.
LABELS = {
    "Measured": r"(?:new\s+)?measurements?",
    "Candidate": r"candidate(?:\s+discover(?:y|ies))?s?",
    "Prediction": r"(?:forecasts?|predictions?)",
    "Correction": r"corrections?",
    "Null result": r"(?:(?:carefully\s+recorded\s+)?searches\s+that\s+found\s+nothing|null\s+results?)",
}
_LABEL = re.compile(r"(?<!\*)\*(" + "|".join(re.escape(name) for name in LABELS) + r")\*(?!\*)")
_BULLET = re.compile(r"^\s*[-*+]\s+")
_TOTAL = re.compile(r"\bhas\s+(\d[\d,]*)\s+entries\b")


def _resolve(base_dir: str, target: str) -> tuple[str, str]:
    """Split a link target into a normalised repository path and its anchor."""
    path, _, anchor = target.partition("#")
    path = unquote(path.split("?", 1)[0])
    if not path:
        return "", unquote(anchor)
    joined = posixpath.normpath(posixpath.join(base_dir, path)) if not path.startswith("/") else posixpath.normpath(path.lstrip("/"))
    return joined, unquote(anchor)


def catalog_entries(tree: RepoTree) -> list[dict]:
    """Labelled bullets of the top-level README: line, label, text and the result folders they cite."""
    if ROOT_README not in tree:
        return []
    prose, _ = markdown.split_code(tree.text(ROOT_README))
    entries = []
    for number, line in enumerate(prose.splitlines(), 1):
        if not _BULLET.match(line):
            continue
        labels = _LABEL.findall(line)
        if not labels:
            continue
        targets = [link.target for link in markdown.links(line) if not link.external]
        folders = []
        for target in targets:
            path, _ = _resolve("", target)
            parts = path.split("/")
            if len(parts) >= 2 and parts[0] == "results" and parts[1] not in folders and tree.has_dir(f"results/{parts[1]}"):
                folders.append(parts[1])
        entries.append({"line": number, "label": labels[-1], "text": line.strip(), "results": folders,
                        "targets": targets})
    return entries


def labels_by_result(tree: RepoTree) -> dict[str, list[str]]:
    """Result folder -> labels of the catalog entries that cite it."""
    table: dict[str, list[str]] = {}
    for entry in catalog_entries(tree):
        for folder in entry["results"]:
            if entry["label"] not in table.setdefault(folder, []):
                table[folder].append(entry["label"])
    return table


def check_links(tree: RepoTree) -> list[Finding]:
    """Every relative Markdown link must resolve inside the tree."""
    findings: list[Finding] = []
    anchor_cache: dict[str, set[str]] = {}
    checked = 0
    for path in tree.paths:
        if not path.lower().endswith(".md"):
            continue
        prose, _ = markdown.split_code(tree.text(path))
        base = posixpath.dirname(path)
        local = 0
        for link in markdown.links(prose):
            if link.external or not link.target:
                continue
            local += 1
            target, anchor = _resolve(base, link.target)
            if target.startswith(".."):
                findings.append(Finding(CHECK, "catalog/link-outside-repository", FAIL, path,
                                        f"link target {link.target!r} leaves the repository", line=link.line,
                                        key=link.target))
                continue
            document = target or path
            if target and target not in tree and not tree.has_dir(target):
                findings.append(Finding(CHECK, "catalog/link-missing", FAIL, path,
                                        f"link target {link.target!r} does not exist (resolved to {target})",
                                        line=link.line, key=link.target, data={"resolved": target}))
                continue
            if anchor and document in tree and document.lower().endswith(".md"):
                if document not in anchor_cache:
                    anchor_cache[document] = markdown.anchors(tree.text(document))
                if anchor.lower() not in anchor_cache[document] and anchor not in anchor_cache[document]:
                    findings.append(Finding(CHECK, "catalog/anchor-missing", WARN, path,
                                            f"no heading or id in {document} produces the anchor #{anchor}",
                                            line=link.line, key=link.target))
        checked += local
    findings.append(Finding(CHECK, "catalog/links-resolve", PASS, ROOT_README if ROOT_README in tree else "",
                            f"{checked} relative Markdown links checked", key="links",
                            data={"links_checked": checked,
                                  "broken": sum(1 for f in findings if f.severity == FAIL)}))
    return findings


def check_coverage(tree: RepoTree) -> list[Finding]:
    """Every result folder should be reachable from the top-level catalog."""
    if ROOT_README not in tree:
        return []
    prose, _ = markdown.split_code(tree.text(ROOT_README))
    cited: set[str] = set()
    for link in markdown.links(prose):
        if link.external:
            continue
        path, _ = _resolve("", link.target)
        parts = path.split("/")
        if len(parts) >= 2 and parts[0] == "results":
            cited.add(parts[1])
    findings = []
    folders = tree.result_dirs()
    for folder in folders:
        if folder not in cited:
            findings.append(Finding(CHECK, "catalog/result-unlisted", WARN, f"results/{folder}/README.md"
                                    if f"results/{folder}/README.md" in tree else f"results/{folder}",
                                    f"results/{folder} is not linked from the top-level README catalog",
                                    key=folder, result=folder))
    findings.append(Finding(CHECK, "catalog/results-listed", PASS, ROOT_README,
                            f"{len(folders) - len(findings)} of {len(folders)} result folders are linked from the catalog",
                            key="coverage", data={"folders": len(folders), "listed": len(folders) - len(findings)}))
    return findings


def check_counts(tree: RepoTree) -> list[Finding]:
    """The summary sentence must state the same counts as the labelled entries."""
    if ROOT_README not in tree:
        return []
    text = tree.text(ROOT_README)
    entries = catalog_entries(tree)
    actual = {name: 0 for name in LABELS}
    for entry in entries:
        actual[entry["label"]] += 1
    findings: list[Finding] = []
    for entry in entries:
        if not entry["targets"]:
            findings.append(Finding(CHECK, "catalog/entry-without-link", WARN, ROOT_README,
                                    "catalog entry carries a label but links to no evidence",
                                    line=entry["line"], key=entry["text"][:80]))
    total = _TOTAL.search(text)
    if not total:
        findings.append(Finding(CHECK, "catalog/no-summary-sentence", INFO, ROOT_README,
                                "no 'has N entries' sentence found, so stated counts were not compared",
                                key="summary", data={"labelled_entries": actual}))
        return findings
    line = text.count("\n", 0, total.start()) + 1
    sentence_end = text.find("\n", total.start())
    sentence = text[total.start():sentence_end if sentence_end != -1 else len(text)]
    stated_total = int(total.group(1).replace(",", ""))
    stated: dict[str, int] = {}
    for name, words in LABELS.items():
        found = re.search(r"(\d[\d,]*)\s+" + words, sentence, flags=re.I)
        if found:
            stated[name] = int(found.group(1).replace(",", ""))
    problems = []
    if stated_total != len(entries):
        problems.append(f"states {stated_total} entries, catalog has {len(entries)}")
    for name, count in stated.items():
        if count != actual[name]:
            problems.append(f"states {count} '{name}', catalog has {actual[name]}")
    data = {"stated_total": stated_total, "stated": stated, "labelled_entries": actual, "entries": len(entries)}
    if problems:
        findings.append(Finding(CHECK, "catalog/count-mismatch", FAIL, ROOT_README,
                                "summary sentence disagrees with the labelled entries: " + "; ".join(problems),
                                line=line, key="summary", data=data))
    else:
        findings.append(Finding(CHECK, "catalog/counts-agree", PASS, ROOT_README,
                                f"summary sentence matches the {len(entries)} labelled entries", line=line,
                                key="summary", data=data))
    return findings


def run(tree: RepoTree) -> list[Finding]:
    return check_links(tree) + check_coverage(tree) + check_counts(tree)
