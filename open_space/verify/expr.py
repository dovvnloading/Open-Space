"""A small, safe expression language for recomputing a claim from committed data.

A claim file says, for example::

    "let": {"t": "table('profiles.csv')"},
    "compute": "pearson(t['stix_50_84'], t['fermi_template'])"

The expression is parsed with Python's own grammar and evaluated by walking the
syntax tree. Only arithmetic, comparisons, indexing, literals, comprehensions and
calls to the functions listed here are allowed. There is no attribute access, no
import and no way to reach the filesystem: ``table`` and ``json`` read committed
bytes from the tree, inside the claim's own result folder.
"""

from __future__ import annotations

import ast
import csv
import gzip
import io
import json as _json
import math
import operator
import posixpath

from . import stats

MAX_STEPS = 20_000_000
MAX_ITEMS = 5_000_000


class ExprError(ValueError):
    """The expression is not allowed, or failed while being evaluated."""


def _number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


class Vec:
    """A column: a sequence with element-wise arithmetic, comparisons and masking."""

    __slots__ = ("data",)

    def __init__(self, data) -> None:
        self.data = list(data)
        if len(self.data) > MAX_ITEMS:
            raise ExprError("vector too large")

    def __len__(self) -> int:
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def __repr__(self) -> str:
        head = ", ".join(repr(item) for item in self.data[:6])
        return f"Vec([{head}{', ...' if len(self.data) > 6 else ''}], n={len(self.data)})"

    def __getitem__(self, key):
        if isinstance(key, Vec):
            if len(key) != len(self):
                raise ExprError(f"mask of length {len(key)} applied to a vector of length {len(self)}")
            return Vec(item for item, keep in zip(self.data, key.data) if keep)
        if isinstance(key, slice):
            return Vec(self.data[key])
        if isinstance(key, bool) or not isinstance(key, int):
            raise ExprError("a vector is indexed by an integer, a slice or a mask")
        return self.data[key]

    def _zip(self, other, function, reflected: bool = False) -> "Vec":
        if isinstance(other, Vec):
            if len(other) != len(self):
                raise ExprError(f"vectors of length {len(self)} and {len(other)} cannot be combined")
            pairs = zip(other.data, self.data) if reflected else zip(self.data, other.data)
            return Vec(function(a, b) for a, b in pairs)
        if reflected:
            return Vec(function(other, item) for item in self.data)
        return Vec(function(item, other) for item in self.data)

    def _unary(self, function) -> "Vec":
        return Vec(function(item) for item in self.data)


def _missing_safe(function):
    """Arithmetic with a missing value (NaN) stays missing; division by zero is missing too."""
    def apply(a, b):
        try:
            return function(a, b)
        except ZeroDivisionError:
            return math.nan
    return apply


_BINARY = {
    ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv,
    ast.FloorDiv: operator.floordiv, ast.Mod: operator.mod,
}
_COMPARE = {
    ast.Eq: operator.eq, ast.NotEq: operator.ne, ast.Lt: operator.lt, ast.LtE: operator.le,
    ast.Gt: operator.gt, ast.GtE: operator.ge,
}


def _power(base, exponent):
    if _number(exponent) and abs(exponent) > 4096:
        raise ExprError("exponent too large")
    if _number(base) and _number(exponent) and base != 0 and abs(exponent) * math.log10(abs(base) + 1e-300) > 400:
        raise ExprError("power would overflow")
    try:
        result = base ** exponent
    except (OverflowError, ZeroDivisionError) as error:
        raise ExprError(f"power failed: {error}") from error
    if isinstance(result, complex):
        raise ExprError("power of a negative number to a fractional exponent")
    return result


def _binary(node_type, left, right):
    for operand in (left, right):
        if isinstance(operand, (str, list, tuple, dict)) and node_type is not ast.Add:
            raise ExprError("text and lists support only + (joining); use vec(...) for arithmetic on a list")
    if node_type is ast.Add and isinstance(left, (str, list, tuple)) and type(left) is not type(right):
        raise ExprError("+ joins two values of the same kind")
    if node_type is ast.Pow:
        function = _power
    elif node_type in (ast.BitAnd, ast.BitOr):
        function = (lambda a, b: bool(a) and bool(b)) if node_type is ast.BitAnd else (lambda a, b: bool(a) or bool(b))
    else:
        function = _missing_safe(_BINARY[node_type])
    if isinstance(left, Vec):
        return left._zip(right, function)
    if isinstance(right, Vec):
        return right._zip(left, function, reflected=True)
    if node_type in (ast.BitAnd, ast.BitOr) and not (isinstance(left, bool) and isinstance(right, bool)):
        raise ExprError("& and | combine masks or truth values")
    return function(left, right)


class Table:
    """A committed table. ``t['column']`` is a Vec; ``t[mask]`` keeps the rows where the mask is true."""

    def __init__(self, names: list[str], columns: dict[str, list], label: str) -> None:
        self.names = names
        self.columns = columns
        self.label = label
        self.rows = len(next(iter(columns.values()))) if columns else 0

    def __len__(self) -> int:
        return self.rows

    def __repr__(self) -> str:
        return f"Table({self.label}, rows={self.rows}, columns={self.names})"

    def __getitem__(self, key):
        if isinstance(key, str):
            if key not in self.columns:
                raise ExprError(f"{self.label} has no column {key!r}; columns are {self.names}")
            return Vec(self.columns[key])
        if isinstance(key, Vec):
            if len(key) != self.rows:
                raise ExprError(f"mask of length {len(key)} applied to {self.rows} rows of {self.label}")
            keep = [index for index, flag in enumerate(key.data) if flag]
            return Table(self.names, {name: [column[index] for index in keep]
                                      for name, column in self.columns.items()}, self.label)
        raise ExprError("a table is indexed by a column name or a row mask")


def _cell(text: str):
    """A table cell as a number when it is one, missing (NaN) when blank, text otherwise."""
    stripped = text.strip()
    if not stripped:
        return math.nan
    try:
        if stripped.lstrip("+-").isdigit():
            return int(stripped)             # whole numbers stay exact: 19-digit catalogue ids do not fit a float
        return float(stripped)
    except ValueError:
        return stripped


def _parse_table(text: str, label: str, delimiter: str | None) -> Table:
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ExprError(f"{label} is empty")
    if delimiter is None:
        delimiter = "\t" if label.lower().split(".gz")[0].endswith((".tsv", ".tab")) else ","
    if delimiter == " ":
        rows = [line.split() for line in lines]
    else:
        rows = list(csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter))
    header = [cell.strip() for cell in rows[0]]
    if all(isinstance(_cell(cell), (int, float)) for cell in header):
        header = [f"column{index + 1}" for index in range(len(header))]
        body = rows
    else:
        body = rows[1:]
    names: list[str] = []
    for index, name in enumerate(header):
        name = name or f"column{index + 1}"
        while name in names:
            name += "_"
        names.append(name)
    columns: dict[str, list] = {name: [] for name in names}
    for row in body:
        for index, name in enumerate(names):
            columns[name].append(_cell(row[index]) if index < len(row) else math.nan)
    for name, column in columns.items():
        if any(isinstance(item, str) for item in column):          # a text column stays text throughout
            columns[name] = ["" if isinstance(item, float) and math.isnan(item) else
                             (item if isinstance(item, str) else _plain(item)) for item in column]
    return Table(names, columns, label)


def _plain(value: float) -> str:
    return str(int(value)) if value == int(value) and abs(value) < 1e15 else repr(value)


def _reduce(function):
    """Wrap a statistic so it accepts a Vec or a list and ignores missing values and text."""
    def apply(values, *rest):
        return function(_floats(values), *rest)
    apply.__name__ = getattr(function, "__name__", "reduce")
    return apply


def _floats(values) -> list[float]:
    if isinstance(values, Vec):
        values = values.data
    if not isinstance(values, (list, tuple)):
        raise ExprError("expected a column or a list of numbers")
    out = []
    for item in values:
        if isinstance(item, bool):
            out.append(1.0 if item else 0.0)
        elif _number(item):
            if math.isfinite(item):
                out.append(float(item))
        else:
            raise ExprError(f"expected numbers, found {item!r}")
    return out


def _pairs(function):
    def apply(x, y, *rest):
        xs = x.data if isinstance(x, Vec) else list(x)
        ys = y.data if isinstance(y, Vec) else list(y)
        return function(xs, ys, *rest)
    apply.__name__ = getattr(function, "__name__", "pairs")
    return apply


def _elementwise(function):
    def apply(value, *rest):
        if isinstance(value, Vec):
            return value._unary(lambda item: function(item, *rest))
        return function(value, *rest)
    apply.__name__ = getattr(function, "__name__", "elementwise")
    return apply


def _guard(function):
    """Mathematical domain errors become missing values instead of aborting a whole column."""
    def apply(value, *rest):
        try:
            return function(value, *rest)
        except (ValueError, ZeroDivisionError, OverflowError):
            return math.nan
    apply.__name__ = getattr(function, "__name__", "guard")
    return apply


def _count(values) -> int:
    """Number of true entries of a mask, or of entries of a column that are not missing."""
    data = values.data if isinstance(values, Vec) else list(values)
    return sum(1 for item in data if item is True or (item is not False and item == item and item != ""))


def _unique(values) -> Vec:
    data = values.data if isinstance(values, Vec) else list(values)
    seen: dict = {}
    for item in data:
        if not (isinstance(item, float) and math.isnan(item)):
            seen.setdefault(item, None)
    return Vec(sorted(seen, key=lambda item: (isinstance(item, str), item)))


def _pointer(document, pointer: str):
    """Resolve a JSON pointer such as '/events/0/delay_s'."""
    node = document
    for part in [piece for piece in pointer.split("/") if piece != ""]:
        part = part.replace("~1", "/").replace("~0", "~")
        if isinstance(node, list):
            node = node[int(part)]
        elif isinstance(node, dict):
            if part not in node:
                raise ExprError(f"JSON pointer {pointer!r}: no key {part!r}")
            node = node[part]
        else:
            raise ExprError(f"JSON pointer {pointer!r} runs past a value")
    return node


def _bounded_range(*arguments) -> list:
    span = range(*arguments)
    if len(span) > MAX_ITEMS:
        raise ExprError("range is too long")
    return list(span)


def _factorial(n: int) -> int:
    if isinstance(n, bool) or not isinstance(n, int) or not 0 <= n <= 5000:
        raise ExprError("factorial takes an integer between 0 and 5000")
    return math.factorial(n)


def _comb(n: int, k: int) -> int:
    for value in (n, k):
        if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= 1_000_000:
            raise ExprError("comb takes integers between 0 and 1,000,000")
    return math.comb(n, k)


def _text_test(method: str):
    """Element-wise text predicate: startswith, endswith or contains."""
    def apply(values, needle):
        if not isinstance(needle, str):
            raise ExprError(f"{method} needs text to look for")

        def one(item) -> bool:
            text = item if isinstance(item, str) else ("" if item != item else str(item))
            if method == "contains":
                return needle in text
            return getattr(text, method)(needle)

        return values._unary(one) if isinstance(values, Vec) else one(values)
    apply.__name__ = method
    return apply


def _vector(*parts) -> Vec:
    if len(parts) == 1 and isinstance(parts[0], (list, tuple, Vec)):
        return Vec(parts[0])
    return Vec(parts)


def _dot(a, b) -> float:
    xs, ys = _floats(a), _floats(b)
    if len(xs) != len(ys):
        raise ExprError("dot product of vectors with different lengths")
    return math.fsum(x * y for x, y in zip(xs, ys))


def _cross(a, b) -> Vec:
    (ax, ay, az), (bx, by, bz) = _floats(a), _floats(b)
    return Vec([ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx])


CONSTANTS = {
    "pi": math.pi, "tau": math.tau, "e": math.e, "nan": math.nan, "inf": math.inf,
    "c": 299792458.0,                     # speed of light, m/s (exact)
    "G": 6.67430e-11,                     # CODATA 2018, m^3 kg^-1 s^-2
    "GM_sun": 1.3271244e20,               # IAU 2015 nominal, m^3 s^-2
    "GM_earth": 3.986004e14,              # IAU 2015 nominal, m^3 s^-2
    "GM_jupiter": 1.2668653e17,           # IAU 2015 nominal, m^3 s^-2
    "au": 149597870700.0,                 # m (exact)
    "pc": 3.0856775814913673e16,          # m
    "ly": 9460730472580800.0,             # m (exact)
    "day": 86400.0, "jyear": 31557600.0,  # s
    "arcsec": math.pi / 648000.0, "mas": math.pi / 648000000.0, "deg": math.pi / 180.0,   # rad
    "k_B": 1.380649e-23, "h_planck": 6.62607015e-34,                                      # SI (exact)
    "R_sun": 6.957e8, "R_earth": 6.3781e6, "R_jupiter": 7.1492e7,                          # IAU nominal, m
    "true": True, "false": False, "none": None,
}

FUNCTIONS = {
    # shape and selection
    "len": len, "count": _count, "unique": _unique, "nunique": lambda values: len(_unique(values)),
    "vec": _vector, "list": lambda values: list(values), "range": _bounded_range,
    "sorted": lambda values: Vec(sorted(_floats(values))), "zip": lambda *columns: [tuple(row) for row in zip(*columns)],
    "first": lambda values: (values.data if isinstance(values, Vec) else list(values))[0],
    "last": lambda values: (values.data if isinstance(values, Vec) else list(values))[-1],
    "any": lambda values: any(values), "all": lambda values: all(values),
    "float": float, "int": int, "str": str, "bool": bool,
    "at": _pointer,
    "startswith": _text_test("startswith"), "endswith": _text_test("endswith"), "contains": _text_test("contains"),
    # reductions
    "sum": _reduce(math.fsum), "mean": _reduce(stats.mean), "median": _reduce(stats.median),
    "min": lambda *args: min(_floats(args[0])) if len(args) == 1 else min(args),
    "max": lambda *args: max(_floats(args[0])) if len(args) == 1 else max(args),
    "std": _reduce(stats.std), "pstd": _reduce(lambda data: stats.std(data, 0)), "var": _reduce(stats.variance),
    "rms": _reduce(stats.rms), "mad": _reduce(stats.mad), "quantile": _reduce(stats.quantile),
    "argmax": lambda values: max(range(len(values)), key=lambda index: values[index]),
    "argmin": lambda values: min(range(len(values)), key=lambda index: values[index]),
    "diff": lambda values: Vec(b - a for a, b in zip(list(values)[:-1], list(values)[1:])),
    "dot": _dot, "cross": _cross, "norm": lambda values: math.sqrt(math.fsum(item * item for item in _floats(values))),
    # element-wise mathematics
    "abs": _elementwise(abs), "sqrt": _elementwise(_guard(math.sqrt)), "exp": _elementwise(_guard(math.exp)),
    "log": _elementwise(_guard(math.log)), "log10": _elementwise(_guard(math.log10)),
    "log2": _elementwise(_guard(math.log2)),
    "sin": _elementwise(math.sin), "cos": _elementwise(math.cos), "tan": _elementwise(math.tan),
    "asin": _elementwise(_guard(math.asin)), "acos": _elementwise(_guard(math.acos)), "atan": _elementwise(math.atan),
    "atan2": math.atan2, "hypot": math.hypot, "radians": _elementwise(math.radians),
    "degrees": _elementwise(math.degrees), "floor": _elementwise(_guard(math.floor)),
    "ceil": _elementwise(_guard(math.ceil)), "round": _elementwise(round), "isnan": _elementwise(lambda item: item != item),
    "erf": _elementwise(math.erf), "erfc": _elementwise(math.erfc), "gamma": _elementwise(_guard(math.gamma)),
    "lgamma": _elementwise(_guard(math.lgamma)),
    # combinatorics and distributions
    "comb": _comb, "factorial": _factorial,
    "norm_sf": stats.norm_sf, "norm_cdf": stats.norm_cdf, "norm_isf": stats.norm_isf,
    "chi2_sf": stats.chi2_sf, "binom_pmf": stats.binom_pmf, "binom_sf": stats.binom_sf, "binom_cdf": stats.binom_cdf,
    "poisson_pmf": stats.poisson_pmf, "poisson_sf": stats.poisson_sf, "poisson_cdf": stats.poisson_cdf,
    # association and fitting
    "pearson": _pairs(stats.pearson), "spearman": _pairs(stats.spearman),
    "slope": _pairs(lambda x, y, weights=None: stats.linear_fit(x, y, weights)[0]),
    "intercept": _pairs(lambda x, y, weights=None: stats.linear_fit(x, y, weights)[1]),
    "wmean": _pairs(stats.weighted_mean),
    "chi2": lambda observed, model, error: stats.chi_square(list(observed), list(model), list(error)),
    # periodicity
    "ls_power": _pairs(stats.lomb_scargle_power), "ls_peak": _pairs(stats.lomb_scargle_peak),
    # sky and time
    "angsep": stats.angular_separation, "unit": lambda ra, dec: Vec(stats.unit_vector(ra, dec)),
    "tai_minus_utc": stats.tai_minus_utc, "tt_minus_utc": stats.tt_minus_utc,
    "jd": stats.julian_date, "mjd": stats.modified_julian_date,
}


class Evaluator:
    """Evaluate expressions against one folder of a tree."""

    def __init__(self, reader=None, base: str = "") -> None:
        self._reader = reader            # callable(path) -> bytes, or None when no data access is offered
        self._base = base
        self._steps = 0
        self.names: dict[str, object] = {}
        self._tables: dict[tuple, Table] = {}
        self._documents: dict[str, object] = {}

    # ---- data access, confined to the claim's folder

    def _bytes(self, path: str) -> tuple[bytes, str]:
        if self._reader is None:
            raise ExprError("this expression has no data folder")
        if not isinstance(path, str) or path.startswith(("/", "\\")) or ":" in path:
            raise ExprError(f"data path {path!r} must be relative to the result folder")
        resolved = posixpath.normpath(posixpath.join(self._base, path.replace("\\", "/")))
        if self._base and not (resolved + "/").startswith(self._base.rstrip("/") + "/"):
            raise ExprError(f"data path {path!r} leaves the result folder")
        try:
            data = self._reader(resolved)
        except KeyError as error:
            raise ExprError(f"no committed file {resolved}") from error
        if resolved.lower().endswith(".gz"):
            data = gzip.decompress(data)
        return data, resolved

    def _table(self, path: str, delimiter: str | None = None) -> Table:
        key = (path, delimiter)
        if key not in self._tables:
            data, resolved = self._bytes(path)
            self._tables[key] = _parse_table(data.decode("utf-8-sig", "replace"), resolved, delimiter)
        return self._tables[key]

    def _json(self, path: str, pointer: str | None = None):
        if path not in self._documents:
            data, _ = self._bytes(path)
            try:
                self._documents[path] = _json.loads(data.decode("utf-8-sig"))
            except ValueError as error:
                raise ExprError(f"{path} is not valid JSON: {error}") from error
        document = self._documents[path]
        return _pointer(document, pointer) if pointer else document

    # ---- evaluation

    def evaluate(self, source: str):
        if not isinstance(source, str) or len(source) > 20_000:
            raise ExprError("expression must be a string of reasonable length")
        try:
            tree = ast.parse(source.strip(), mode="eval")
        except SyntaxError as error:
            raise ExprError(f"cannot parse expression: {error.msg}") from error
        try:
            return self._eval(tree.body, {})
        except ExprError:
            raise
        except RecursionError as error:
            raise ExprError("expression nests too deeply") from error
        except Exception as error:          # an error inside a library function is the claim's failure, not ours
            raise ExprError(f"{type(error).__name__}: {error}") from error

    def define(self, name: str, source: str) -> None:
        if not name.isidentifier() or name in FUNCTIONS or name in CONSTANTS or name in ("table", "json"):
            raise ExprError(f"{name!r} cannot be used as a name")
        self.names[name] = self.evaluate(source)

    def _tick(self) -> None:
        self._steps += 1
        if self._steps > MAX_STEPS:
            raise ExprError("expression is too expensive to evaluate")

    def _eval(self, node, scope: dict):
        self._tick()
        if isinstance(node, ast.Constant):
            if isinstance(node.value, (int, float, str, bool)) or node.value is None:
                return node.value
            raise ExprError(f"literal {node.value!r} is not allowed")
        if isinstance(node, ast.Name):
            if node.id in scope:
                return scope[node.id]
            if node.id in self.names:
                return self.names[node.id]
            if node.id in CONSTANTS:
                return CONSTANTS[node.id]
            raise ExprError(f"unknown name {node.id!r}")
        if isinstance(node, ast.BinOp):
            kind = type(node.op)
            if kind not in _BINARY and kind not in (ast.Pow, ast.BitAnd, ast.BitOr):
                raise ExprError(f"operator {kind.__name__} is not allowed")
            return _binary(kind, self._eval(node.left, scope), self._eval(node.right, scope))
        if isinstance(node, ast.UnaryOp):
            value = self._eval(node.operand, scope)
            if isinstance(node.op, ast.USub):
                return value._unary(operator.neg) if isinstance(value, Vec) else -value
            if isinstance(node.op, ast.UAdd):
                return value
            if isinstance(node.op, (ast.Not, ast.Invert)):
                return value._unary(lambda item: not item) if isinstance(value, Vec) else (not value)
            raise ExprError("unary operator is not allowed")
        if isinstance(node, ast.BoolOp):
            result = None
            for part in node.values:
                result = self._eval(part, scope)
                if isinstance(result, Vec):
                    raise ExprError("use & and | to combine masks; 'and' and 'or' need single truth values")
                if isinstance(node.op, ast.And) and not result:
                    return result
                if isinstance(node.op, ast.Or) and result:
                    return result
            return result
        if isinstance(node, ast.Compare):
            left = self._eval(node.left, scope)
            outcome = None
            for op, comparator in zip(node.ops, node.comparators):
                right = self._eval(comparator, scope)
                if type(op) not in _COMPARE:
                    raise ExprError(f"comparison {type(op).__name__} is not allowed")
                function = _COMPARE[type(op)]
                if isinstance(left, Vec):
                    step = left._zip(right, function)
                elif isinstance(right, Vec):
                    step = right._zip(left, function, reflected=True)
                else:
                    step = function(left, right)
                if outcome is None:
                    outcome = step
                elif isinstance(outcome, Vec) or isinstance(step, Vec):
                    outcome = _binary(ast.BitAnd, outcome, step)
                else:
                    outcome = outcome and step
                left = right
            return outcome
        if isinstance(node, ast.IfExp):
            return self._eval(node.body if self._eval(node.test, scope) else node.orelse, scope)
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name):
                raise ExprError("only the listed functions can be called")
            if any(isinstance(argument, ast.Starred) for argument in node.args) or any(
                    keyword.arg is None for keyword in node.keywords):
                raise ExprError("argument unpacking is not allowed")
            arguments = [self._eval(argument, scope) for argument in node.args]
            keywords = {keyword.arg: self._eval(keyword.value, scope) for keyword in node.keywords}
            name = node.func.id
            if name == "table":
                return self._table(*arguments, **keywords)
            if name == "json":
                return self._json(*arguments, **keywords)
            if name not in FUNCTIONS:
                raise ExprError(f"unknown function {name!r}")
            return FUNCTIONS[name](*arguments, **keywords)
        if isinstance(node, ast.Subscript):
            container = self._eval(node.value, scope)
            if isinstance(node.slice, ast.Slice):
                bounds = [None if part is None else self._eval(part, scope)
                          for part in (node.slice.lower, node.slice.upper, node.slice.step)]
                key = slice(*bounds)
            else:
                key = self._eval(node.slice, scope)
            if isinstance(container, (Table, Vec)):
                return container[key]
            if isinstance(container, (list, tuple, str)):
                if isinstance(key, bool) or not isinstance(key, (int, slice)):
                    raise ExprError("sequences are indexed by integers or slices")
                return container[key]
            if isinstance(container, dict):
                if key not in container:
                    raise ExprError(f"no key {key!r}")
                return container[key]
            raise ExprError(f"{type(container).__name__} cannot be indexed")
        if isinstance(node, (ast.List, ast.Tuple)):
            items = [self._eval(item, scope) for item in node.elts]
            return items if isinstance(node, ast.List) else tuple(items)
        if isinstance(node, ast.Dict):
            if any(key is None for key in node.keys):
                raise ExprError("dictionary unpacking is not allowed")
            return {self._eval(key, scope): self._eval(value, scope) for key, value in zip(node.keys, node.values)}
        if isinstance(node, (ast.ListComp, ast.GeneratorExp)):
            results: list = []
            self._comprehend(node.generators, 0, node.elt, scope, results)
            return results
        raise ExprError(f"{type(node).__name__} is not allowed in a claim expression")

    def _comprehend(self, generators, depth: int, element, scope: dict, results: list) -> None:
        generator = generators[depth]
        if generator.is_async:
            raise ExprError("async comprehensions are not allowed")
        iterable = self._eval(generator.iter, scope)
        if isinstance(iterable, Table):
            raise ExprError("iterate over a column, not a whole table")
        if isinstance(iterable, dict):
            iterable = list(iterable)
        for item in iterable:
            self._tick()
            inner = dict(scope)
            self._bind(generator.target, item, inner)
            if all(self._eval(condition, inner) for condition in generator.ifs):
                if depth + 1 < len(generators):
                    self._comprehend(generators, depth + 1, element, inner, results)
                else:
                    results.append(self._eval(element, inner))
                    if len(results) > MAX_ITEMS:
                        raise ExprError("comprehension produces too many items")

    def _bind(self, target, value, scope: dict) -> None:
        if isinstance(target, ast.Name):
            scope[target.id] = value
        elif isinstance(target, (ast.Tuple, ast.List)):
            values = list(value)
            if len(values) != len(target.elts):
                raise ExprError("cannot unpack: wrong number of values")
            for part, item in zip(target.elts, values):
                self._bind(part, item, scope)
        else:
            raise ExprError("comprehension target must be a name or a tuple of names")


def evaluate(source: str, names: dict | None = None):
    """Evaluate a pure expression with no data access. Convenience for tests and identities."""
    evaluator = Evaluator()
    evaluator.names.update(names or {})
    return evaluator.evaluate(source)
