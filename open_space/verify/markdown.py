"""Small Markdown reader: code regions, links, headings and line numbers.

It is not a full parser. It recognises exactly the constructs the checks need,
and it keeps every character offset stable so that findings can cite a line.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass
import re

_FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
_INLINE_CODE = re.compile(r"(?<!`)(`+)(?!`)([^\n]*?)(?<!`)\1(?!`)")
_LINK = re.compile(r"(!?)\[((?:[^\[\]]|\[[^\[\]]*\])*)\]\(\s*(<[^<>\n]*>|[^()\s]*(?:\([^()\s]*\)[^()\s]*)*)(?:\s+(?:\"[^\"\n]*\"|'[^'\n]*'))?\s*\)")
_HEADING = re.compile(r"^ {0,3}(#{1,6})[ \t]+(.*?)[ \t]*#*[ \t]*$")
_HTML_ID = re.compile(r"<[a-zA-Z][^<>]*?\b(?:id|name)\s*=\s*[\"']([^\"']+)[\"']")
_SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.\-]*:")


class Lines:
    """Translate character offsets into 1-based line numbers."""

    def __init__(self, text: str) -> None:
        self._starts = [0]
        for found in re.finditer("\n", text):
            self._starts.append(found.end())

    def at(self, offset: int) -> int:
        return bisect.bisect_right(self._starts, offset)


def _blank(chars: list[str], start: int, end: int) -> None:
    for index in range(start, end):
        if chars[index] != "\n":
            chars[index] = " "


@dataclass
class CodeBlock:
    line: int          # line of the first content line
    language: str
    content: str
    fenced: bool


def split_code(text: str) -> tuple[str, list[CodeBlock]]:
    """Return the prose with code blanked out (same length), plus the code that was removed."""
    chars = list(text)
    blocks: list[CodeBlock] = []
    offset = 0
    open_fence: tuple[str, int, str, int, int] | None = None     # char, length, language, start offset, line
    line_number = 0
    for raw in text.splitlines(keepends=True):
        line_number += 1
        stripped = raw.rstrip("\r\n")
        match = _FENCE.match(stripped)
        if open_fence is None:
            if match:
                marker = match.group(1)
                open_fence = (marker[0], len(marker), match.group(2).strip().split(" ")[0].lower(),
                              offset + len(raw), line_number + 1)
                _blank(chars, offset, offset + len(raw))
        else:
            char, length, language, start, first_line = open_fence
            closing = match and match.group(1)[0] == char and len(match.group(1)) >= length and not match.group(2).strip()
            if closing:
                blocks.append(CodeBlock(first_line, language, text[start:offset], True))
                open_fence = None
            _blank(chars, offset, offset + len(raw))
        offset += len(raw)
    if open_fence is not None:                                   # unterminated fence runs to the end
        char, length, language, start, first_line = open_fence
        blocks.append(CodeBlock(first_line, language, text[start:], True))
    prose = "".join(chars)
    lines = Lines(text)
    chars = list(prose)
    for match in _INLINE_CODE.finditer(prose):
        blocks.append(CodeBlock(lines.at(match.start()), "", match.group(2).strip(), False))
        _blank(chars, match.start(), match.end())
    blocks.sort(key=lambda block: (block.line, not block.fenced))
    return "".join(chars), blocks


@dataclass
class Link:
    line: int
    text: str
    target: str        # as written, without angle brackets
    image: bool
    start: int         # character span of the whole construct in the document
    end: int

    @property
    def external(self) -> bool:
        return bool(_SCHEME.match(self.target)) or self.target.startswith("//")


def links(prose: str) -> list[Link]:
    """Inline links and images in text whose code has already been blanked."""
    lines = Lines(prose)
    found = []
    for match in _LINK.finditer(prose):
        target = match.group(3)
        if target.startswith("<") and target.endswith(">"):
            target = target[1:-1]
        found.append(Link(lines.at(match.start()), match.group(2), target, bool(match.group(1)),
                          match.start(), match.end()))
    return found


def slug(heading: str) -> str:
    """GitHub's anchor for a heading: lower case, punctuation dropped, spaces to hyphens."""
    text = re.sub(r"!\[([^\]]*)\]\([^)]*\)", r"\1", heading)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"<[^>]+>", "", text)
    kept = [char for char in text.strip().lower() if char.isalnum() or char in "-_ "]
    return "".join(kept).replace(" ", "-")


def anchors(text: str) -> set[str]:
    """Every anchor a document offers: heading slugs (with duplicates numbered) and explicit ids."""
    prose, _ = split_code(text)
    offered: set[str] = set()
    counts: dict[str, int] = {}
    for line in prose.splitlines():
        match = _HEADING.match(line)
        if not match:
            continue
        for base in {slug(match.group(2)), slug(match.group(2).replace("_", "").replace("*", ""))}:
            seen = counts.get(base, 0)
            counts[base] = seen + 1
            offered.add(base if seen == 0 else f"{base}-{seen}")
    offered.update(_HTML_ID.findall(text))
    return offered
