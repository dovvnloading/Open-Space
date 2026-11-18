"""Build the visitor catalog from committed evidence. Run from any directory."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
import subprocess
import sys
import tarfile
from urllib.parse import quote, unquote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from open_space.verify import markdown
from open_space.verify.catalog import catalog_entries, check_counts
from open_space.verify.model import FAIL
from open_space.verify.runner import DEFAULT_CHECKS, LEDGER_SCHEMA
from open_space.verify.tree import RepoTree
from observatory.research_destinations import build_destinations
from observatory.verify_spotlights import verify_spotlights

# These two shared-report entries were originally published with title hashes.
# Preserve their public URLs across revised titles and future catalog entries.
PUBLISHED_SHARED_IDS = {
    ("results/neowise-nuclear-flares/README.md", ""): "results-neowise-nuclear-flares-28374e29",
    ("results/neowise-nuclear-flares/README.md", "no-vanishing-supergiant"):
        "results-neowise-nuclear-flares-3dc81494",
}


def git(*arguments: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(ROOT), *arguments])


def read_reports(revision: str) -> RepoTree:
    """Index all tracked paths, reading only the small report documents as blobs."""
    names = git("ls-tree", "-r", "-z", "--name-only", revision).decode().split("\0")
    files = {name: b"" for name in names if name}
    packed = git("archive", "--format=tar", revision, "--", ":(glob)**/*.md")
    with tarfile.open(fileobj=io.BytesIO(packed)) as archive:
        for item in archive:
            if item.isfile():
                files[item.name] = archive.extractfile(item).read()
    return RepoTree(files, commit=revision)


def publication_times(revision: str) -> dict[str, str]:
    """Date when a report change appeared on the revision's first-parent history."""
    log = git("log", "--first-parent", "-m", "--no-renames", "--name-only",
              "--format=@%H%x09%cI", revision, "--", "results", "docs").decode()
    shallow_file = Path(git("rev-parse", "--path-format=absolute", "--git-path", "shallow").decode().strip())
    boundaries = set(shallow_file.read_text().split()) if shallow_file.exists() else set()
    dates = {}
    timestamp = None
    for line in log.splitlines():
        if line.startswith("@"):
            commit, timestamp = line[1:].split("\t")
            if commit in boundaries:
                timestamp = None
        elif timestamp and line.startswith(("results/", "docs/")):
            dates.setdefault(line, timestamp)
    return dates


def plain(text: str) -> str:
    for link in reversed(markdown.links(text)):
        text = text[:link.start] + link.text + text[link.end:]
    return " ".join(text.replace("**", "").replace("`", "").strip().split())


def report_heading(tree: RepoTree, path: str) -> str:
    if path.endswith(".md") and path in tree:
        for line in tree.text(path).splitlines():
            if line.startswith("# "):
                return plain(line[2:])
    return PurePosixPath(path).name.replace("-", " ").replace("_", " ")


def report_introduction(tree: RepoTree, path: str) -> str:
    """Retain the source's complete opening paragraph, including qualifications."""
    if path not in tree or not path.endswith(".md"):
        return ""
    paragraphs = re.split(r"\n\s*\n", tree.text(path).strip())
    for paragraph in paragraphs:
        if paragraph.startswith(("#", "!", "<!--")):
            continue
        if paragraph.startswith(("|", "- ", "* ", "```")):
            return ""  # An index or table does not have a narrative abstract.
        return plain(paragraph)
    return ""


def source_link(tree: RepoTree, repository: str, target: str, label: str) -> dict:
    path, _, anchor = target.partition("#")
    path = posixpath.normpath(unquote(path).lstrip("/") or "README.md")
    anchor = unquote(anchor)
    kind = "blob" if path in tree else "tree"
    if kind == "tree" and not tree.has_dir(path):
        raise ValueError(f"Missing evidence at {tree.commit}: {target}")
    url = f"{repository}/{kind}/{tree.commit}/{quote(path)}"
    if anchor:
        document = path if kind == "blob" else f"{path}/README.md"
        if document not in tree or not document.lower().endswith(".md"):
            raise ValueError(f"Unsupported evidence anchor at {tree.commit}: {target}")
        if anchor not in markdown.anchors(tree.text(document)):
            raise ValueError(f"Missing evidence anchor at {tree.commit}: {target}")
        url += "#" + quote(anchor)
    return {"label": plain(label), "path": path, "url": url}


def record_id(path: str) -> str:
    parts = list(PurePosixPath(path).with_suffix("").parts)
    if parts[-1].lower() == "readme":
        parts.pop()
    return "-".join(parts).lower().replace("_", "-")


def build_records(tree: RepoTree, repository: str, dates: dict[str, str]) -> list[dict]:
    for finding in check_counts(tree):
        if finding.severity == FAIL:
            raise ValueError(f"Source catalog disagrees with its entries: {finding.message}")
    records = []
    covered_folders = set()
    lines = tree.text("README.md").splitlines()
    topics = {}
    topic = "Research archive"
    for number, line in enumerate(lines, 1):
        if line.startswith("### "):
            topic = plain(line[4:])
        topics[number] = topic
    for entry in catalog_entries(tree):
        line = lines[entry["line"] - 1]
        match = re.match(r"^\s*[-*+]\s+\*\*(.+?)\*\*\s*(.*)", line)
        if not match:
            raise ValueError(f"Expected a titled catalog entry on README line {entry['line']}")
        title, description = match.groups()
        links = [source_link(tree, repository, link.target, link.text)
                 for link in markdown.links(line) if not link.external and not link.image]
        if not links:
            raise ValueError(f"No local evidence on README line {entry['line']}")
        path = links[0]["path"]
        summary = plain(description.split(f"*{entry['label']}*", 1)[0])
        records.append({
            "title": plain(title).removesuffix("."), "summary": summary or report_heading(tree, path),
            "status": entry["label"], "topic": topics[entry["line"]], "links": links,
            "provenance": {"kind": "catalog",
                           "url": f"{repository}/blob/{tree.commit}/README.md?plain=1#L{entry['line']}"},
        })
        covered_folders.update(entry["results"])

    for folder in tree.result_dirs():
        if folder in covered_folders:
            continue
        directory = f"results/{folder}"
        readme = directory + "/README.md"
        documents = [p for p in tree.under(directory) if p.endswith(".md") and p.count("/") == 2]
        path = readme if readme in tree else next(iter(documents), directory)
        link = source_link(tree, repository, path, "Research record")
        records.append({
            "title": report_heading(tree, path), "summary": report_introduction(tree, path),
            "status": "Research record",
            "topic": "Research archive", "links": [link],
            "provenance": {"kind": "report", "url": link["url"]},
        })

    ids = Counter(record_id(record["links"][0]["path"]) for record in records)
    for record in records:
        main = record["links"][0]
        key = PUBLISHED_SHARED_IDS.get((main["path"], unquote(main["url"].partition("#")[2])))
        if key is None:
            key = record_id(main["path"])
            if ids[key] > 1:
                key += "-" + hashlib.sha256(record["title"].encode()).hexdigest()[:8]
        updates = [dates[link["path"]] for link in record["links"] if link["path"] in dates]
        record.update(id=key, path=main["path"], url=main["url"],
                      updatedAt=max(updates, key=datetime.fromisoformat, default=None))
    if len({r["id"] for r in records}) != len(records):
        raise ValueError("Duplicate observatory record IDs")
    return records


def archive_audit(tree: RepoTree, repository: str, records: list[dict], ledger: dict | None) -> dict:
    """Attach scoped verifier results only when the ledger describes this committed tree."""
    links = [link for row in records for link in row["links"]]
    classified = sum(row["provenance"]["kind"] == "catalog" for row in records)
    audit = {"catalogEntries": classified, "supplementalRecords": len(records) - classified,
             "recordLinks": len(links), "uniqueEvidencePaths": len({link["path"] for link in links}),
             "repositoryChecks": None}
    for row in records:
        row["evidenceChecks"] = []
    if ledger is None:
        return audit

    identity = ledger["tree"]
    if (ledger["schema"] != LEDGER_SCHEMA or identity["origin"] != "commit"
            or identity["commit"] != tree.commit or ledger["scope"] != "all"
            or not set(DEFAULT_CHECKS).issubset(ledger["checks"])):
        raise ValueError("Expected a whole-repository verification ledger for this committed revision")
    summary = ledger["summary"]
    if summary["findings"]["fail"]:
        raise ValueError("Verification ledger contains failed checks")
    cards = {card["result"]: card for card in ledger["results"]}
    if set(cards) != set(tree.result_dirs()) or len(cards) != len(ledger["results"]):
        raise ValueError("Verification ledger does not cover each result folder exactly once")

    checksums, claims, numbers = summary["checksums"], summary["claims"], summary["numbers"]
    catalog_trace = numbers["catalog"]
    audit["repositoryChecks"] = {
        "matchedDigests": checksums["verified_against_committed_bytes"],
        "externalReceipts": checksums["external_receipts"], "claimFiles": claims["claim_files"],
        "declaredClaims": claims["claims"], "reproducedClaims": claims["reproduced"],
        "reportNumbers": numbers["precise"], "tracedReportNumbers": numbers["found"],
        "untracedReportNumbers": numbers["unmatched"], "catalogTraceEntries": catalog_trace["entries"],
        "preciseCatalogNumbers": catalog_trace["precise"],
        "tracedCatalogNumbers": catalog_trace["confirmed"] + catalog_trace["found"],
        "openWarnings": summary["findings"]["warn"], "resultFolders": len(cards),
        "methodUrl": source_link(tree, repository, "docs/VERIFICATION.md", "Verification methods")["url"],
    }
    for row in records:
        # A grouped finding can cite several folders. Attribute the verifier's
        # counts to those folders without transferring them to the finding itself.
        locations = {link["path"].split("/")[1] for link in row["links"] if link["path"].startswith("results/")}
        for folder in sorted(locations.intersection(cards)):
            card = cards[folder]
            numbers = card["numbers"] or {}
            row["evidenceChecks"].append({
                "folder": folder, "url": source_link(tree, repository, f"results/{folder}", folder)["url"],
                "matchedDigests": card["checksums"]["verified"], "externalReceipts": card["checksums"]["external"],
                "declaredClaims": card["claims"]["declared"], "reproducedClaims": card["claims"]["reproduced"],
                "preciseNumbers": numbers.get("precise", 0), "tracedNumbers": numbers.get("found", 0),
                "openWarnings": card["warn"],
            })
    return audit


def build_spotlights(tree: RepoTree, repository: str, records: list[dict], output: Path) -> list[dict]:
    curated = json.loads(git("show", f"{tree.commit}:observatory/data/spotlights.json"))
    # Older published manifests predate typed readouts. Republish their original
    # text without assigning a numerical replay that was never performed.
    typed = any("key" in metric for item in curated for metric in item["measurements"])
    verify_readouts = "observatory/verify_spotlights.py" in tree or typed
    if verify_readouts:
        verify_spotlights(curated, lambda path: git("show", f"{tree.commit}:{path}"))
    for item in curated:
        path = item["path"]
        matches = [r for r in records if any(link["path"] == path for link in r["links"])]
        if len(matches) != 1:
            raise ValueError(f"Expected a single catalog entry for spotlight {path}")
        item["recordId"] = matches[0]["id"]
        item["status"] = matches[0]["status"]
        item["updatedAt"] = matches[0]["updatedAt"]
        main = source_link(tree, repository, path, "Investigation report")
        item["url"] = main["url"]
        item["links"] = [main] + [source_link(tree, repository, p, label) for label, p in item.pop("evidence")]
        if verify_readouts:
            for measurement in item["measurements"]:
                measurement["sources"] = [source_link(tree, repository, p, PurePosixPath(p).name)
                                          for p in measurement.pop("sourcePaths")]
            table_path = item["table"].pop("sourcePath")
            item["table"]["source"] = source_link(tree, repository, table_path, PurePosixPath(table_path).name)
        figure = item["figure"]
        figure["url"] = source_link(tree, repository, figure["path"], "Figure source")["url"]
        name = PurePosixPath(figure["src"]).name
        destination = output.parent / "figures" / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if figure["path"].endswith((".svg", ".png")):
            destination.write_bytes(git("show", f"{tree.commit}:{figure['path']}"))
        else:
            image = git("show", f"{tree.commit}:observatory/data/figures/{name}")
            # Matplotlib embeds the digest of the exact CSV drawn. Prevent an
            # archived figure from silently illustrating a changed data table.
            digest = hashlib.sha256(git("show", f"{tree.commit}:{figure['path']}")).hexdigest()
            if digest.encode() not in image:
                raise ValueError(f"Regenerate {name} from {figure['path']} with render_figures.py")
            source = HERE / "data/figures" / name
            if destination.resolve() == source.resolve() and source.exists():
                if source.read_bytes() != image:
                    raise ValueError(f"Use --output outside the source directory to preserve {name}")
            else:
                destination.write_bytes(image)
    return curated


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", default="HEAD", help="Committed Git revision to publish")
    parser.add_argument("--repository", default="dovvnloading/Open-Space", help="GitHub owner/repository")
    parser.add_argument("--output", type=Path, default=HERE / "data/catalog.json")
    parser.add_argument("--verification-ledger", type=Path,
                        help="Include folder-level checks from a full open_space.verify JSON ledger at this revision")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repository):
        parser.error("--repository must be a GitHub owner/repository")
    revision = git("rev-parse", "--verify", "--end-of-options", args.revision + "^{commit}").decode().strip()
    repository = f"https://github.com/{args.repository}"
    tree = read_reports(revision)
    records = build_records(tree, repository, publication_times(revision))
    ledger = json.loads(args.verification_ledger.read_bytes()) if args.verification_ledger else None
    audit = archive_audit(tree, repository, records, ledger)
    spotlights = build_spotlights(tree, repository, records, args.output)
    destinations = build_destinations(tree, records, spotlights, repository, args.output,
                                      lambda path: git("show", f"{tree.commit}:{path}"))
    catalog = {
        "schemaVersion": 1, "revision": revision, "repository": repository,
        "revisionDate": git("show", "-s", "--format=%cI", revision).decode().strip(),
        "records": records, "spotlights": spotlights, "destinations": destinations, "audit": audit,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"Published {len(records)} research records and {len(spotlights)} observing dossiers from {revision[:10]}")


if __name__ == "__main__":
    main()
