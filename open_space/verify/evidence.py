"""The machine-readable values deposited in a result folder.

Four pools, searched in this order:

* direct: JSON values and the cells of small tables. This is where a result
  states its numbers.
* derived: facts computed from the files rather than stored in them: row counts,
  counts per category, extremes, sums, means, medians, list lengths, and numbers
  embedded in short text.
* bulk: the cells of large tables and long text files, the raw rows.
* code: numeric constants in the folder's scripts, the parameters of the analysis.

Each value remembers where it was first seen, so a match can be inspected by hand.
"""

from __future__ import annotations

import bisect
from collections import Counter
import csv
import gzip
import io
from itertools import compress, zip_longest
import json
import math
import re
import statistics

from .numbers import decimals_of
from .tree import RepoTree

SMALL_TABLE_ROWS = 50          # tables up to this size state results; larger ones hold raw rows
SMALL_DOCUMENT_LEAVES = 3000
SHORT_TEXT_LINES = 400
MAX_FILE_BYTES = 60_000_000
LOW_CARDINALITY = 12
CROSS_TAB_COLUMNS = 6
UNKNOWN = 99                   # decimal places of a value that was computed rather than written
CLAIMS_FILE = "claims.json"    # declares what the prose says; reading it as evidence would be circular

TABLE_SUFFIXES = (".csv", ".tsv", ".tab")
TEXT_SUFFIXES = (".txt", ".dat", ".tbl", ".par", ".mpc", ".lc", ".out", ".log", ".ecsv", ".list", ".lis")
_FLOAT = re.compile(r"(?<![\w.])[+\-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+\-]?\d+)?(?![\w]|\.\d)")
_NOISE = re.compile(r"\d{4}-\d{2}-\d{2}(?:[T ]\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?)?Z?"
                    r"|(?<![0-9A-Za-z])[0-9a-fA-F]{24,}(?![0-9A-Za-z])")
_POSITIVE = (0.0).__lt__
_NEGATIVE = (0.0).__gt__

# A lexer just large enough to find numeric literals in Python source: comments, strings and
# names are recognised so that digits inside them are skipped. Python's own parser is not used,
# because a script written for a newer Python would then be read on some machines and not others.
_PYTHON = re.compile(r"""
    (?P<comment>\#[^\n]*)
  | (?P<string>[rRbBuUfF]{0,3}(?:\"\"\"(?:\\.|[^\\])*?\"\"\"|'''(?:\\.|[^\\])*?'''
                               |\"(?:\\.|[^\"\\\n])*\"|'(?:\\.|[^'\\\n])*'))
  | (?P<name>[A-Za-z_][A-Za-z_0-9]*)
  | (?P<number>(?:0[xX][0-9a-fA-F_]+|0[oO][0-7_]+|0[bB][01_]+
                |(?:\d[\d_]*(?:\.[\d_]*)?|\.\d[\d_]*)(?:[eE][+\-]?\d[\d_]*)?)[jJ]?)
""", re.VERBOSE | re.DOTALL)
_BEFORE_UNARY = set("([{,=:+-*/%<>!&|^~;")
_KEYWORDS_BEFORE_UNARY = {"return", "in", "and", "or", "not", "if", "else", "elif", "yield", "lambda", "is", "while",
                          "assert", "print"}


def code_numbers(text: str):
    """Yield (value, line) for each real-valued numeric literal in Python source, with its sign."""
    line = 1
    cursor = 0
    previous = ""                       # last significant token: an operator character, a name or "number"
    for match in _PYTHON.finditer(text):
        gap = text[cursor:match.start()]
        line += gap.count("\n")
        stripped = gap.strip()
        minus = False
        if stripped:
            if stripped.endswith("-"):
                head = stripped[:-1].rstrip()
                # A minus is a sign when an operator, an opening bracket or a keyword comes before it.
                minus = head[-1] in _BEFORE_UNARY if head else (previous in _KEYWORDS_BEFORE_UNARY or previous == "")
            previous = stripped[-1]
        kind = match.lastgroup
        if kind == "number":
            token = match.group().replace("_", "")
            if token[-1] not in "jJ":
                try:
                    value = float(int(token, 0)) if token[:2].lower() in ("0x", "0o", "0b") else float(token)
                except (ValueError, OverflowError):
                    value = None
                if value is not None and value == value and value not in (math.inf, -math.inf):
                    yield (-value if minus else value), line
            previous = "number"
        elif kind == "name":
            previous = match.group()
        elif kind == "string":
            previous = "number"         # a string closes an operand, like a number does
        line += match.group().count("\n")
        cursor = match.end()


class _Literal(float):
    """A float that remembers how many decimal places its text had."""

    __slots__ = ("places",)


def _literal(text: str) -> _Literal:
    number = _Literal(text)
    number.places = decimals_of(text)
    return number


def _to_float(cell: str) -> float | None:
    if not cell or len(cell) > 40:
        return None
    try:
        value = float(cell)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _plain_places(cell: str) -> int:
    """Decimal places of a cell already known to be numeric. Exponent forms take the slow path."""
    if "e" in cell or "E" in cell:
        places = decimals_of(cell)
        return UNKNOWN if places is None else places
    dot = cell.find(".")
    return 0 if dot < 0 else len(cell) - dot - 1


def _tokens(text: str):
    """Numbers embedded in free text, as (value, decimal places)."""
    for token in _FLOAT.findall(_NOISE.sub(" ", text)):
        try:
            value = float(token.replace("d", "e").replace("D", "e"))
        except ValueError:
            continue
        if math.isfinite(value):
            places = decimals_of(token)
            yield value, UNKNOWN if places is None else places


class Pool:
    """Sorted distinct values. Each one keeps the place it was first seen.

    Values arrive one at a time (``add``) or as a whole table column (``add_column``).
    Columns are merged with dictionary operations only, so millions of raw cells cost
    little; where a column value came from is worked out when someone asks.
    """

    def __init__(self, name: str) -> None:
        self.name = name
        self._first: dict[float, tuple] = {}     # value -> (places, source index, detail)
        self._columns: list[tuple[dict, int, str]] = []
        self._rows: dict[float, int] = {}
        self._owner: dict[float, int] = {}
        self.values: list[float] = []

    def add(self, value: float, source: int, detail, places: int | None = None) -> None:
        places = UNKNOWN if places is None else places
        known = self._first.get(value)
        if known is None:
            self._first[value] = (places, source, detail)
        elif places < known[0]:                  # remember the least precise spelling of this value
            self._first[value] = (places,) + known[1:]

    def add_column(self, first_row: dict[float, int], source: int, name: str) -> None:
        """Register a column as {value: first row}. Values must be finite."""
        self._columns.append((first_row, source, name))

    def freeze(self) -> None:
        for index in range(len(self._columns) - 1, -1, -1):      # earlier columns win
            column = self._columns[index][0]
            self._rows.update(column)
            self._owner.update(dict.fromkeys(column, index))
        self.values = sorted(self._first.keys() | self._rows.keys())

    def __len__(self) -> int:
        return len(self.values)

    def _meta(self, index: int) -> tuple:
        value = self.values[index]
        known = self._first.get(value)
        if known is not None:
            return known
        _, source, name = self._columns[self._owner[value]]
        return (UNKNOWN, source, (name, self._rows[value]))

    def places(self, index: int) -> int | None:
        """Decimal places the value was written with, or None when it was computed or unknown."""
        places = self._meta(index)[0]
        return None if places == UNKNOWN else places

    def source(self, index: int) -> int:
        return self._meta(index)[1]

    def detail(self, index: int):
        return self._meta(index)[2]

    def nearest(self, target: float) -> int | None:
        if not self.values:
            return None
        at = bisect.bisect_left(self.values, target)
        best = None
        for index in (at - 1, at):
            if 0 <= index < len(self.values):
                if best is None or abs(self.values[index] - target) < abs(self.values[best] - target):
                    best = index
        return best


class Evidence:
    """Everything numeric under one folder, plus its raw text for exact identifier lookups."""

    def __init__(self, folder: str) -> None:
        self.folder = folder
        self.sources: list[str] = []
        self.direct = Pool("direct")
        self.derived = Pool("derived")
        self.bulk = Pool("bulk")
        self.code = Pool("code")
        self.texts: list[tuple[str, str]] = []
        self.tables = 0
        self.documents = 0
        self.scripts = 0
        self.skipped: list[str] = []

    @property
    def pools(self) -> tuple[Pool, Pool, Pool, Pool]:
        return self.direct, self.derived, self.bulk, self.code

    def describe(self, pool: Pool, index: int) -> str:
        """Human-readable location of a value: file plus pointer, column and row, or line."""
        source = self.sources[pool.source(index)]
        detail = pool.detail(index)
        if isinstance(detail, tuple):
            column, row = detail
            return f"{source} [{column}] row {row}"
        return f"{source} {detail}".rstrip()

    def find_literal(self, needle: str) -> str | None:
        """Path of the first evidence file whose text contains ``needle`` exactly."""
        for path, text in self.texts:
            if needle in text:
                return path
        return None

    def _source(self, relative: str) -> int:
        self.sources.append(relative)
        return len(self.sources) - 1

    # ------------------------------------------------------------------ builders

    @classmethod
    def build(cls, tree: RepoTree, folder: str) -> "Evidence":
        evidence = cls(folder)
        prefix = folder.rstrip("/") + "/"
        for path in tree.under(folder):
            relative = path[len(prefix):]
            lowered = path.lower()
            if lowered.endswith((".md", ".markdown")) or lowered.endswith("/" + CLAIMS_FILE) or lowered == CLAIMS_FILE:
                continue                                  # prose and claims are what is being checked, not evidence
            compressed = lowered.endswith(".gz")
            if compressed:
                lowered = lowered[:-3]
            if not lowered.endswith(TABLE_SUFFIXES + TEXT_SUFFIXES + (".json", ".geojson", ".py")):
                continue
            data = tree.read(path)
            if compressed:
                try:
                    data = gzip.decompress(data)
                except (OSError, EOFError, ValueError):
                    evidence.skipped.append(relative)
                    continue
            if len(data) > MAX_FILE_BYTES or b"\0" in data[:4096]:
                evidence.skipped.append(relative)
                continue
            text = data.decode("utf-8-sig", "replace")
            if lowered.endswith((".json", ".geojson")):
                evidence._add_json(relative, text)
            elif lowered.endswith(TABLE_SUFFIXES):
                evidence._add_table(relative, text, "\t" if lowered.endswith((".tsv", ".tab")) else ",")
            elif lowered.endswith(".py"):
                evidence._add_code(relative, text)
            else:
                evidence._add_text(relative, text)
            evidence.texts.append((relative, text))
        for pool in evidence.pools:
            pool.freeze()
        return evidence

    def _add_json(self, relative: str, text: str) -> None:
        try:
            document = json.loads(text, parse_float=_literal, parse_int=_literal)
        except (ValueError, RecursionError):
            self._add_text(relative, text)
            return
        self.documents += 1
        source = self._source(relative)
        leaves: list[tuple[float, int | None, str]] = []
        stack = [(document, "")]
        while stack:
            node, pointer = stack.pop()
            if isinstance(node, dict):
                if len(node) > 1:
                    self.derived.add(float(len(node)), source, f"#{pointer or '/'} (number of keys)", 0)
                for key, value in node.items():
                    stack.append((value, f"{pointer}/{key}"))
            elif isinstance(node, list):
                self.derived.add(float(len(node)), source, f"#{pointer or '/'} (number of items)", 0)
                for index, value in enumerate(node):
                    stack.append((value, f"{pointer}/{index}"))
            elif isinstance(node, bool) or node is None:
                continue
            elif isinstance(node, float):
                if math.isfinite(node):
                    leaves.append((float(node), getattr(node, "places", None), pointer))
            elif isinstance(node, str):
                whole = _to_float(node.strip())
                if whole is not None:
                    leaves.append((whole, decimals_of(node), pointer))
                elif len(node) <= 2000:
                    for count, (value, places) in enumerate(_tokens(node)):
                        if count >= 200:
                            break
                        self.derived.add(value, source, f"#{pointer} (inside text)", places)
        pool = self.direct if len(leaves) <= SMALL_DOCUMENT_LEAVES else self.bulk
        for value, places, pointer in leaves:
            pool.add(value, source, "#" + (pointer or "/"), places)

    def _add_table(self, relative: str, text: str, delimiter: str) -> None:
        try:
            rows = [row for row in csv.reader(io.StringIO(text), delimiter=delimiter)
                    if row and not row[0].lstrip().startswith("#")]
        except csv.Error:
            self._add_text(relative, text)
            return
        if len(rows) < 2 or max(map(len, rows)) < 2:
            self._add_text(relative, text)
            return
        self.tables += 1
        source = self._source(relative)
        header = [cell.strip() for cell in rows[0]]
        if all(_to_float(cell) is not None for cell in header if cell):
            body = rows                                   # no header row: every row is data
            header = [f"column {index + 1}" for index in range(len(header))]
        else:
            body = rows[1:]
        derived = self.derived
        count = len(body)
        derived.add(float(count), source, "(number of rows)", 0)
        small = count <= SMALL_TABLE_ROWS
        columns = list(zip_longest(*body, fillvalue=""))[:len(header)]
        categorical: list[tuple[str, tuple]] = []
        for column, name in enumerate(header):
            if column >= len(columns):
                break
            name = name or f"column {column + 1}"
            cells = columns[column]
            blanks = cells.count("")
            numeric: list[float] = []
            row_numbers = None
            try:                                           # a fully numeric column parses in one C-level pass
                if blanks:
                    numeric = list(map(float, filter(None, cells)))
                    row_numbers = list(compress(range(1, count + 1), cells))
                else:
                    numeric = list(map(float, cells))
                    row_numbers = range(1, count + 1)
            except ValueError:
                cells = tuple(cell.strip() for cell in cells)
                blanks = cells.count("")
                numeric, row_numbers = [], []
                for row_number, cell in enumerate(cells, 1):
                    value = _to_float(cell)
                    if value is not None:
                        numeric.append(value)
                        row_numbers.append(row_number)
            filled = count - blanks
            distinct = None
            if numeric:
                first_row = dict(zip(reversed(numeric), reversed(row_numbers)))
                if not math.isfinite(math.fsum(map(abs, first_row))):          # NaN or infinity present
                    first_row = {value: row for value, row in first_row.items() if math.isfinite(value)}
                    numeric = [value for value in numeric if math.isfinite(value)]
                if small:
                    for value, row_number in first_row.items():
                        self.direct.add(value, source, (name, row_number), _plain_places(cells[row_number - 1].strip()))
                else:
                    self.bulk.add_column(first_row, source, name)
                if len(numeric) == filled:
                    distinct = len(first_row)
            if filled and blanks:
                derived.add(float(filled), source, f"[{name}] (non-empty cells)", 0)
            if distinct is None or distinct <= LOW_CARDINALITY:
                labels = Counter(cells)
                labels.pop("", None)
                distinct = len(labels)
                if 1 < distinct <= LOW_CARDINALITY and filled > distinct:
                    for label, occurrences in labels.items():
                        derived.add(float(occurrences), source, f"[{name}] (rows equal to {label!r})", 0)
                    if len(categorical) < CROSS_TAB_COLUMNS:
                        categorical.append((name, cells))
            if 1 < distinct < filled:
                derived.add(float(distinct), source, f"[{name}] (distinct values)", 0)
            if len(numeric) >= 2 and len(numeric) >= 0.8 * filled:
                total = math.fsum(numeric)
                derived.add(min(numeric), source, f"[{name}] (minimum)")
                derived.add(max(numeric), source, f"[{name}] (maximum)")
                derived.add(total, source, f"[{name}] (sum)")
                derived.add(total / len(numeric), source, f"[{name}] (mean)")
                derived.add(statistics.median(numeric), source, f"[{name}] (median)")
                positive = sum(map(_POSITIVE, numeric))
                negative = sum(map(_NEGATIVE, numeric))
                if positive and negative:
                    derived.add(float(positive), source, f"[{name}] (positive values)", 0)
                    derived.add(float(negative), source, f"[{name}] (negative values)", 0)
        # Counts for pairs of categories: "610 g-band rows pass the flag cut".
        for first in range(len(categorical)):
            for second in range(first + 1, len(categorical)):
                name_a, cells_a = categorical[first]
                name_b, cells_b = categorical[second]
                for (a, b), occurrences in Counter(zip(cells_a, cells_b)).items():
                    if a and b:
                        derived.add(float(occurrences), source, f"(rows with {name_a}={a!r} and {name_b}={b!r})", 0)

    def _add_text(self, relative: str, text: str) -> None:
        source = self._source(relative)
        pool = self.derived if text.count("\n") < SHORT_TEXT_LINES else self.bulk
        for number, line in enumerate(text.splitlines(), 1):
            for value, places in _tokens(line):
                pool.add(value, source, f"line {number}", places)

    def _add_code(self, relative: str, text: str) -> None:
        """Numeric literals of a script. The script is neither run nor compiled."""
        self.scripts += 1
        source = self._source(relative)
        for number, line in code_numbers(text):
            self.code.add(number, source, f"line {line} (constant)")
