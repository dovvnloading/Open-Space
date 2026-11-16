"""Checksum checks: every SHA-256 recorded for committed bytes must equal those bytes.

A recorded digest is found in four kinds of place: checksum manifests
(``SHA256SUMS`` style), JSON documents, tables with a digest column, and free
text. Each one is resolved to the file it speaks about when the surrounding
structure names one, and to any committed file with those bytes otherwise.
Digests of bytes that were never committed (upstream downloads) are counted as
external receipts: they are recorded, and they cannot be checked offline.

A mismatch is never reported bare. The check first tries to explain it: line
ending conversion, a byte-order mark, a stale digest that matches an earlier
commit of the file, or a manifest entry that describes a different, larger
upstream object with the same name.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import posixpath
import re
import zipfile

from .model import FAIL, INFO, PASS, WARN, Finding, result_of
from .tree import RepoTree

CHECK = "integrity"
MAX_TEXT_BYTES = 12_000_000

_HEX = re.compile(r"(?<![0-9A-Za-z])[0-9a-fA-F]{64}(?![0-9A-Za-z])")
_HEX_FULL = re.compile(r"[0-9a-fA-F]{64}\Z")
_MANIFEST_LINE = re.compile(r"^([0-9a-fA-F]{64})[ \t]+\*?(\S.*?)\s*$")
_URL = re.compile(r"^(?:https?|ftp)://", re.I)

PATH_KEYS = ("path", "file", "filename", "file_name", "local", "local_path", "relative_path", "relpath",
             "output", "artifact", "saved_as", "name", "target", "deposit", "product")
SIZE_KEYS = ("bytes", "size", "size_bytes", "nbytes", "length", "content_length", "filesize", "file_size")
URL_KEYS = ("url", "source_url", "uri", "href", "download_url", "source", "link")

# A digest key says which bytes it describes. "sha256" next to a file name means that file.
# "raw_sha256", "archive_sha256" and the like describe an upstream object, not the committed file.
_DIGEST_SUFFIX = re.compile(r"[_\-]?(?:sha[_\-]?256(?:sum|_hex|_digest)?|digest|checksum|hash)$")
_PATH_SUFFIX = re.compile(r"[_\-](?:file|path|name|filename)$")


def _digest_stem(key: str) -> str | None:
    """'' for a plain digest key, the qualifier for 'csv_sha256' style keys, None if not a digest key."""
    lowered = key.lower()
    if not _DIGEST_SUFFIX.search(lowered):
        return None
    return _DIGEST_SUFFIX.sub("", lowered)


def _bind(stem: str | None, siblings: list[tuple[str, str]], parent: str | None) -> str | None:
    """Choose the committed file a digest key speaks about, or None when it names no local file."""
    if stem is None:
        return None
    if stem == "":
        # Only a key that names the hashed object counts. A sibling such as "used_for" or
        # "extraction" may mention a committed file without the digest being that file's.
        for want in PATH_KEYS:
            for name, resolved in siblings:
                if name == want:
                    return resolved
        return parent
    for name, resolved in siblings:
        if name == stem or _PATH_SUFFIX.sub("", name) == stem:
            return resolved
    return None


def _is_text(data: bytes) -> bool:
    return len(data) <= MAX_TEXT_BYTES and b"\0" not in data[:8192]


class _Resolver:
    """Turn a path written inside some file into a path of the tree, if it exists there."""

    def __init__(self, tree: RepoTree) -> None:
        self.tree = tree

    def resolve(self, source: str, written: str) -> str | None:
        if not isinstance(written, str) or not written or len(written) > 400 or "\n" in written:
            return None
        if _URL.match(written):
            return None
        written = written.replace("\\", "/")
        if written.startswith("./"):
            written = written[2:]
        bases = [posixpath.dirname(source)]
        owner = result_of(source)
        if owner:
            bases.append(f"results/{owner}")
        bases.append("")
        for base in bases:
            candidate = posixpath.normpath(posixpath.join(base, written))
            if candidate.startswith("..") or candidate == ".":
                continue
            if candidate in self.tree:
                return candidate
        return None


class _Assertion:
    """One recorded digest and what is known about the bytes it describes."""

    __slots__ = ("source", "line", "digest", "written", "bound", "url", "size", "kind")

    def __init__(self, source, line, digest, written=None, bound=None, url=None, size=None, kind="text"):
        self.source = source
        self.line = line
        self.digest = digest.lower()
        self.written = written      # path as written next to the digest
        self.bound = bound          # tree path the structure ties the digest to
        self.url = url
        self.size = size
        self.kind = kind            # manifest, json, table or text


def _line_of(text: str, needle: str, start: int = 0) -> int | None:
    at = text.find(needle, start)
    return None if at < 0 else text.count("\n", 0, at) + 1


def _first(mapping: dict, keys: tuple[str, ...], kind) -> object | None:
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for key in keys:
        value = lowered.get(key)
        if isinstance(value, kind) and not isinstance(value, bool):
            return value
    return None


def _manifest_assertions(path: str, text: str, resolver: _Resolver) -> list[_Assertion] | None:
    """Parse a ``<digest>  <path>`` listing. Returns None when the file is not such a listing."""
    rows = []
    meaningful = 0
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        meaningful += 1
        match = _MANIFEST_LINE.match(line)
        if match:
            rows.append((number, match.group(1), match.group(2)))
    if not rows or len(rows) < 0.8 * meaningful:
        return None
    return [_Assertion(path, number, digest, written=written, bound=resolver.resolve(path, written), kind="manifest")
            for number, digest, written in rows]


def _json_assertions(path: str, text: str, document, resolver: _Resolver) -> list[_Assertion]:
    found: list[_Assertion] = []
    cursor = {"at": 0}

    def locate(digest: str) -> int | None:
        at = text.find(digest, cursor["at"])
        if at < 0:
            at = text.find(digest)
        if at < 0:
            return None
        cursor["at"] = at + 1
        return text.count("\n", 0, at) + 1

    def walk(node, parent_key: str | None) -> None:
        if isinstance(node, dict):
            url = _first(node, URL_KEYS, str)
            if url is not None and not _URL.match(url):
                url = None
            size = _first(node, SIZE_KEYS, int)
            siblings: list[tuple[str, str]] = []            # (key, resolved path) for path-like values
            for key, value in node.items():
                if isinstance(value, str) and not _HEX_FULL.match(value):
                    resolved = resolver.resolve(path, value)
                    if resolved:
                        siblings.append((str(key).lower(), resolved))
            parent_path = resolver.resolve(path, parent_key) if parent_key is not None else None
            for key, value in node.items():
                if isinstance(value, str) and _HEX_FULL.match(value):
                    key_text = str(key)
                    direct = resolver.resolve(path, key_text)
                    if direct:                               # {"data/x.csv": "<digest>"}
                        bound = direct
                    else:                                    # {"file": "x.csv", "sha256": ...} or {"x.csv": {"sha256": ...}}
                        bound = _bind(_digest_stem(key_text), siblings, parent_path)
                    plain = direct is not None or _digest_stem(key_text) == ""
                    found.append(_Assertion(path, locate(value), value, written=bound, bound=bound, url=url,
                                            size=size if plain else None, kind="json"))
                else:
                    walk(value, str(key))
        elif isinstance(node, list):
            for value in node:
                walk(value, parent_key)
        elif isinstance(node, str):
            for match in _HEX.finditer(node):                # digest embedded in a longer string
                if not _HEX_FULL.match(node):
                    found.append(_Assertion(path, locate(match.group()), match.group(), kind="json"))

    walk(document, None)
    return found


def _table_assertions(path: str, text: str, resolver: _Resolver) -> list[_Assertion]:
    delimiter = "\t" if path.lower().endswith((".tsv", ".tab")) else ","
    found: list[_Assertion] = []
    try:
        rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    except csv.Error:
        return _text_assertions(path, text)
    if not rows:
        return found
    header = [cell.strip().lower() for cell in rows[0]]
    for number, row in enumerate(rows, 1):
        digests = [(index, cell.strip()) for index, cell in enumerate(row) if _HEX_FULL.match(cell.strip())]
        if not digests:
            continue
        url = next((cell.strip() for cell in row if _URL.match(cell.strip())), None)
        siblings = []
        for index, cell in enumerate(row):
            if index < len(header) and not _HEX_FULL.match(cell.strip()):
                resolved = resolver.resolve(path, cell.strip())
                if resolved:
                    siblings.append((header[index], resolved))
        line = _line_of(text, digests[0][1]) or number
        for index, digest in digests:
            # The column name decides which bytes a digest describes, exactly as a JSON key does.
            stem = _digest_stem(header[index]) if index < len(header) else None
            bound = _bind(stem, siblings, None)
            found.append(_Assertion(path, line, digest, written=bound, bound=bound, url=url, kind="table"))
    return found


def _text_assertions(path: str, text: str) -> list[_Assertion]:
    return [_Assertion(path, text.count("\n", 0, match.start()) + 1, match.group(), kind="text")
            for match in _HEX.finditer(text)]


def collect(tree: RepoTree) -> list[_Assertion]:
    """Every recorded SHA-256 in the tree, tied to a committed file where the structure allows."""
    resolver = _Resolver(tree)
    assertions: list[_Assertion] = []
    for path in tree.paths:
        data = tree.read(path)
        if not _is_text(data) or not _HEX.search(data.decode("ascii", "replace")):
            continue
        text = tree.text(path)
        lowered = path.lower()
        manifest = _manifest_assertions(path, text, resolver)
        if manifest is not None:
            assertions.extend(manifest)
        elif lowered.endswith(".json"):
            try:
                document = json.loads(text)
            except ValueError:
                assertions.extend(_text_assertions(path, text))
            else:
                assertions.extend(_json_assertions(path, text, document, resolver))
        elif lowered.endswith((".csv", ".tsv", ".tab")):
            assertions.extend(_table_assertions(path, text, resolver))
        else:
            assertions.extend(_text_assertions(path, text))
    return assertions


def content_index(tree: RepoTree) -> dict[str, str]:
    """Digest -> where those bytes live: a file, a ZIP member or a gzip payload."""
    index: dict[str, str] = {}
    for digest, paths in tree.by_digest().items():
        index[digest] = paths[0]
    for path in tree.paths:
        lowered = path.lower()
        if lowered.endswith(".zip"):
            try:
                with zipfile.ZipFile(io.BytesIO(tree.read(path))) as bundle:
                    for member in bundle.infolist():
                        if not member.is_dir() and member.file_size <= 200_000_000:
                            index.setdefault(hashlib.sha256(bundle.read(member)).hexdigest(),
                                             f"{path}!{member.filename}")
            except (zipfile.BadZipFile, OSError, RuntimeError):
                continue
        elif lowered.endswith(".gz"):
            try:
                payload = gzip.decompress(tree.read(path))
            except (OSError, EOFError, ValueError):
                continue
            index.setdefault(hashlib.sha256(payload).hexdigest(), f"{path} (decompressed)")
    return index


def _variants(data: bytes):
    """Byte transformations that commonly explain a digest computed on a different copy."""
    unix = data.replace(b"\r\n", b"\n")
    yield "CRLF line endings (a Windows checkout)", unix.replace(b"\n", b"\r\n")
    yield "LF line endings", unix
    yield "no final newline", data.rstrip(b"\r\n")
    yield "one added final newline", data.rstrip(b"\r\n") + b"\n"
    if data.startswith(b"\xef\xbb\xbf"):
        yield "no byte-order mark", data[3:]
    else:
        yield "a byte-order mark", b"\xef\xbb\xbf" + data
    if data[:2] == b"\x1f\x8b":
        try:
            yield "its decompressed payload", gzip.decompress(data)
        except (OSError, EOFError, ValueError):
            pass


def explain(tree: RepoTree, assertion: _Assertion) -> tuple[str, str, str, dict]:
    """Classify a mismatch. Returns (rule, severity, message, data)."""
    path = assertion.bound
    data = tree.read(path)
    actual = tree.digest(path)
    facts = {"recorded": assertion.digest, "committed": actual, "file": path, "committed_bytes": len(data)}
    if assertion.size is not None:
        facts["recorded_bytes"] = assertion.size
    for label, variant in _variants(data):
        if variant != data and hashlib.sha256(variant).hexdigest() == assertion.digest:
            facts["explained_by"] = label
            return ("integrity/digest-of-converted-bytes", WARN,
                    f"recorded digest of {path} equals the committed file with {label}, not the committed bytes", facts)
    for commit, digest in tree.history(path):
        if digest == assertion.digest and digest != actual:
            facts["matches_commit"] = commit
            return ("integrity/stale-digest", FAIL,
                    f"recorded digest of {path} matches the file as of commit {commit[:10]}; "
                    "the file changed afterwards and the digest was not updated", facts)
    if assertion.size is not None and assertion.size != len(data):
        if assertion.url:
            facts["url"] = assertion.url
            return ("integrity/name-shadows-upstream-file", WARN,
                    f"the entry for {assertion.written!r} describes an upstream object of {assertion.size:,} bytes; "
                    f"the committed {path} is a different file of {len(data):,} bytes with no digest of its own",
                    facts)
        return ("integrity/digest-mismatch", FAIL,
                f"recorded digest and size ({assertion.size:,} bytes) do not match committed {path} "
                f"({len(data):,} bytes)", facts)
    return ("integrity/digest-mismatch", FAIL,
            f"recorded digest {assertion.digest[:16]}… does not match committed {path} ({actual[:16]}…)", facts)


def run(tree: RepoTree) -> tuple[list[Finding], dict]:
    """Verify every recorded digest. Returns findings and a census for the report."""
    assertions = collect(tree)
    index = content_index(tree)
    findings: list[Finding] = []
    per_source: dict[str, dict] = {}
    external: dict[str, dict] = {}
    census = {"assertions": len(assertions), "distinct_digests": len({a.digest for a in assertions}),
              "path_bound_verified": 0, "content_verified": 0, "external": 0, "listed_but_absent": 0,
              "mismatched": 0}

    def tally(source: str) -> dict:
        return per_source.setdefault(source, {"verified": 0, "content": 0, "external": 0, "absent": 0})

    # Files whose committed bytes carry a digest of their own somewhere in the repository.
    vouched = {assertion.bound for assertion in assertions
               if assertion.bound and tree.digest(assertion.bound) == assertion.digest}

    for assertion in assertions:
        counts = tally(assertion.source)
        if assertion.bound:
            if tree.digest(assertion.bound) == assertion.digest:
                counts["verified"] += 1
                census["path_bound_verified"] += 1
                continue
            if assertion.digest in index and assertion.kind != "manifest":
                # The digest belongs to other committed bytes; the path association was a guess.
                counts["content"] += 1
                census["content_verified"] += 1
                continue
            rule, severity, message, data = explain(tree, assertion)
            if rule != "integrity/name-shadows-upstream-file" or assertion.bound not in vouched:
                census["mismatched"] += 1
                findings.append(Finding(CHECK, rule, severity, assertion.source, message, line=assertion.line,
                                        key=f"{assertion.digest}:{assertion.bound}", data=data))
                continue
            # An upstream download that shares its name with a committed file. The committed file has
            # its own verified digest, so the two cannot be confused: this one is an external receipt.
        if assertion.digest in index:
            counts["content"] += 1
            census["content_verified"] += 1
            holder = index[assertion.digest]
            inside_archive = "!" in holder and assertion.written and \
                holder.split("!", 1)[1].endswith(assertion.written.replace("\\", "/").lstrip("./"))
            if assertion.kind == "manifest" and assertion.written and not inside_archive:
                findings.append(Finding(CHECK, "integrity/manifest-path-moved", WARN, assertion.source,
                                        f"manifest lists {assertion.written!r}, which is not in the repository; "
                                        f"identical bytes are at {index[assertion.digest]}", line=assertion.line,
                                        key=f"{assertion.digest}:{assertion.written}",
                                        data={"listed": assertion.written, "found_at": index[assertion.digest]}))
            continue
        if assertion.kind == "manifest":
            counts["absent"] += 1
            census["listed_but_absent"] += 1
        else:
            counts["external"] += 1
        record = external.setdefault(assertion.digest, {"sources": []})
        if assertion.source not in record["sources"]:
            record["sources"].append(assertion.source)
        if assertion.url and "url" not in record:
            record["url"] = assertion.url
        if assertion.size is not None and "bytes" not in record:
            record["bytes"] = assertion.size

    census["external"] = len(external)
    census["external_with_url"] = sum(1 for record in external.values() if "url" in record)
    for source in sorted(per_source):
        counts = per_source[source]
        checked = counts["verified"] + counts["content"]
        if checked:
            findings.append(Finding(CHECK, "integrity/digests-verified", PASS, source,
                                    f"{checked} recorded digests equal committed bytes", key="verified",
                                    data={"path_bound": counts["verified"], "by_content": counts["content"]}))
        if counts["absent"]:
            findings.append(Finding(CHECK, "integrity/listed-file-not-deposited", INFO, source,
                                    f"{counts['absent']} manifest entries name files that are not in the repository",
                                    key="absent", data={"entries": counts["absent"]}))
        if counts["external"]:
            findings.append(Finding(CHECK, "integrity/external-receipts", INFO, source,
                                    f"{counts['external']} digests describe bytes outside the repository "
                                    "(recorded, not checkable offline)", key="external",
                                    data={"digests": counts["external"]}))
    return findings, {"census": census, "external": external}
