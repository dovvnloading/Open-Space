"""Build the visitor catalog from committed evidence. Run from any directory."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tarfile
from urllib.parse import quote, unquote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from open_space.verify import markdown
from open_space.verify.catalog import catalog_entries
from open_space.verify.tree import RepoTree


def git(*arguments: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(ROOT), *arguments])


def read_reports(revision: str) -> RepoTree:
    """Index all tracked paths, reading only the small report documents as blobs."""
    names = git("ls-tree", "-r", "-z", "--name-only", revision).decode().split("\0")
    files = {name: b"" for name in names if name}
    packed = git("archive", "--format=tar", revision, "--", "README.md",
                 ":(glob)docs/*.md", ":(glob)results/**/*.md")
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
    path = unquote(path).removeprefix("./").rstrip("/")
    kind = "blob" if path in tree else "tree"
    if kind == "tree" and not tree.has_dir(path):
        raise ValueError(f"Missing evidence at {tree.commit}: {target}")
    url = f"{repository}/{kind}/{tree.commit}/{quote(path)}"
    if anchor:
        url += "#" + quote(unquote(anchor))
    return {"label": plain(label), "path": path, "url": url}


def record_id(path: str) -> str:
    parts = list(PurePosixPath(path).with_suffix("").parts)
    if parts[-1].lower() == "readme":
        parts.pop()
    return "-".join(parts).lower().replace("_", "-")


def build_records(tree: RepoTree, repository: str, dates: dict[str, str]) -> list[dict]:
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
        })
        covered_folders.update(entry["results"])

    for folder in tree.result_dirs():
        if folder in covered_folders:
            continue
        directory = f"results/{folder}"
        readme = directory + "/README.md"
        documents = [p for p in tree.under(directory) if p.endswith(".md") and p.count("/") == 2]
        path = readme if readme in tree else next(iter(documents), directory)
        records.append({
            "title": report_heading(tree, path), "summary": report_introduction(tree, path),
            "status": "Research record",
            "topic": "Research archive", "links": [source_link(tree, repository, path, "Research record")],
        })

    ids = Counter(record_id(record["links"][0]["path"]) for record in records)
    for record in records:
        main = record["links"][0]
        key = record_id(main["path"])
        if ids[key] > 1:
            key += "-" + hashlib.sha256(record["title"].encode()).hexdigest()[:8]
        updates = [dates[link["path"]] for link in record["links"] if link["path"] in dates]
        record.update(id=key, path=main["path"], url=main["url"],
                      updatedAt=max(updates, key=datetime.fromisoformat, default=None))
    if len({r["id"] for r in records}) != len(records):
        raise ValueError("Duplicate observatory record IDs")
    return records


def build_spotlights(tree: RepoTree, repository: str, records: list[dict], output: Path) -> list[dict]:
    curated = json.loads((HERE / "data/spotlights.json").read_text(encoding="utf-8"))
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
        figure = item["figure"]
        figure["url"] = source_link(tree, repository, figure["path"], "Figure source")["url"]
        name = PurePosixPath(figure["src"]).name
        destination = output.parent / "figures" / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if figure["path"].endswith((".svg", ".png")):
            destination.write_bytes(git("show", f"{tree.commit}:{figure['path']}"))
        else:
            image = (HERE / "data/figures" / name).read_bytes()
            # Matplotlib embeds the digest of the exact CSV drawn. Prevent an
            # archived figure from silently illustrating a changed data table.
            digest = hashlib.sha256(git("show", f"{tree.commit}:{figure['path']}")).hexdigest()
            if digest.encode() not in image:
                raise ValueError(f"Regenerate {name} from {figure['path']} with render_figures.py")
            if destination.resolve() != (HERE / "data/figures" / name).resolve():
                destination.write_bytes(image)
    return curated


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", default="HEAD", help="Committed Git revision to publish")
    parser.add_argument("--repository", default="dovvnloading/Open-Space", help="GitHub owner/repository")
    parser.add_argument("--output", type=Path, default=HERE / "data/catalog.json")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repository):
        parser.error("--repository must be a GitHub owner/repository")
    revision = git("rev-parse", "--verify", "--end-of-options", args.revision + "^{commit}").decode().strip()
    repository = f"https://github.com/{args.repository}"
    tree = read_reports(revision)
    records = build_records(tree, repository, publication_times(revision))
    spotlights = build_spotlights(tree, repository, records, args.output)
    catalog = {
        "schemaVersion": 1, "revision": revision, "repository": repository,
        "revisionDate": git("show", "-s", "--format=%cI", revision).decode().strip(),
        "records": records, "spotlights": spotlights,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"Published {len(records)} research records and {len(spotlights)} observing dossiers from {revision[:10]}")


if __name__ == "__main__":
    main()
