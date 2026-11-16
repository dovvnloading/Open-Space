"""Read the numbers a Markdown document states, with the precision each one claims.

A printed number is a claim about an interval, not a point: "12.99" asserts a
value in [12.985, 12.995]. This module finds every number in the prose, works
out that interval, and separates measurements from the things that merely look
like numbers: catalogue designations, years, dates, table and section labels.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import re

from . import markdown

SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
MINUS = "−"

_NUMBER = re.compile(r"""
    (?P<sign>[+\-−])?
    (?P<int>\d{1,3}(?:,\d{3})+|\d+)
    (?P<frac>\.\d+)?
    (?:
        \s?(?:×|\\times|\\cdot|·|[xX*])\s?10
        (?:\^\{?\s?(?P<e1>[+\-−]?\d+)\s?\}?|(?P<e2>[⁺⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))
      |
        [eE](?P<e3>[+\-−]?\d+)(?![A-Za-z0-9_.])
    )?
""", re.VERBOSE)

_POWER_OF_TEN = re.compile(r"\^\{?\s?([+\-−]?\d+)\s?\}?|([⁺⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)")
_TEN_TO = re.compile(r"(?:×|\\times|\\cdot|·)?\s?10(?:\^\{?\s?([+\-−]?\d+)\s?\}?|([⁺⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))")

_URL = re.compile(r"(?:https?|ftp)://[^\s<>)\]]+")
_HTML = re.compile(r"<[^<>\n]{1,200}>")
_HEX_RUN = re.compile(r"(?<![0-9A-Za-z])(?=[0-9a-fA-F]*[a-fA-F])(?=[a-fA-F]*[0-9])[0-9a-fA-F]{12,}(?![0-9A-Za-z])")
_ISO_DATE = re.compile(r"(?<![\d.])\d{4}-\d{2}-\d{2}(?:[T ]\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?)?(?:Z|\s?UTC)?")
_SLASH_DATE = re.compile(r"(?<![\d.])\d{4}/\d{1,2}/\d{1,2}(?![\d])")
_CLOCK = re.compile(r"(?<![\d.:])\d{1,2}:\d{2}(?::\d{2}(?:\.\d+)?)?(?![\d:])")
_SEXAGESIMAL = re.compile(r"(?<![\w.])[+\-−]?\d{1,3}[hd°]\s?\d{1,2}[m′']\s?\d{1,2}(?:\.\d+)?[s″\"]?")

_MONTH = r"(?:Jan(?:uary)?|Feb(?:ruary)?|Mar(?:ch)?|Apr(?:il)?|May|Jun(?:e)?|Jul(?:y)?|Aug(?:ust)?|Sep(?:t(?:ember)?)?|Oct(?:ober)?|Nov(?:ember)?|Dec(?:ember)?)"
_MONTH_BEFORE = re.compile(r"\b" + _MONTH + r"\.?,?\s$")
_MONTH_AFTER = re.compile(r"^(?:st|nd|rd|th)?,?\s" + _MONTH + r"\b")
_DAY_RANGE_AFTER = re.compile(r"^\s?(?:[–\-]|and|,|to)\s?\d{1,2}(?:\s?(?:[–\-]|and|,|to)\s?\d{1,2})*,?\s" + _MONTH + r"\b")

_REFERENCE = re.compile(
    r"(?:\b(?:Tables?|Figures?|Figs?\.?|Equations?|Eqs?\.?|Sections?|Sects?\.?|Appendix|Chapters?|Refs?\.?|"
    r"versions?|Release|Data\s+Release|DR|EDR|PR|Issue|issue|step|Step|Note|Notes|No\.|Vol\.|vol\.|pp?\.)|[§#])\s?$")
# A catalogue prefix directly before an integer: "TIC 384244835", "HD 34445", "2dFS 757".
# Time scales and common measurement acronyms are quantities, not catalogues.
_DESIGNATION = re.compile(
    r"(?<![A-Za-z0-9])(?!(?:MJD|JD|BJD|HJD|BTJD|TJD|TBJD|TDB|UTC|TT|TAI|RA|DEC|SNR|RMS|FWHM|AU|PA|ID|OF|AND|TO|"
    r"IS|IN|AT|BY|ON|OR|THE)[ \-]?$)"
    r"(?:(?=[A-Z0-9]{2,}[ \-]?$)[A-Z0-9]*[A-Z][A-Z0-9]*|[0-9]+[A-Za-z]+[A-Z][A-Za-z]*)[ \-]?$")
_APPROX = re.compile(r"(?:\b(?:about|approximately|approx\.?|around|roughly|nearly|almost|some|circa|ca\.)|[~≈∼≃])\s?$",
                     re.I)
_MULTIPLIER = re.compile(r"^\s(thousand|million|billion|trillion)\b", re.I)
_MULTIPLIERS = {"thousand": 3, "million": 6, "billion": 9, "trillion": 12}
_CONNECTOR = re.compile(r"^\s*(?:and|to|or|through|–|—|-|,|±|\\pm|and\s+about|to\s+about)\s*[+\-−]?$")
_TABLE_ROW = re.compile(r"^\s*\|")
_TABLE_RULE = re.compile(r"^\s*\|?\s*:?-{2,}:?\s*(?:\|\s*:?-{2,}:?\s*)*\|?\s*$")

GLUE = "/+-" + MINUS + "–.:"


@dataclass
class Number:
    """One number printed in a document."""

    text: str                  # exactly as printed
    value: float
    ulp: float                 # size of one unit in the last printed place
    sig: int                   # significant digits claimed
    kind: str                  # value, identifier, year, date or reference
    line: int
    start: int
    end: int
    percent: bool = False
    approximate: bool = False
    signed: bool = False
    context: str = ""
    readings: list[tuple[float, float, str]] = field(default_factory=list)   # alternative (value, ulp, why)

    @property
    def half(self) -> float:
        return self.ulp / 2

    @property
    def precise(self) -> bool:
        """At least three significant digits: a wrong number rarely matches by accident."""
        return self.kind == "value" and self.sig >= 3

    def interval(self) -> tuple[float, float]:
        return self.value - self.half, self.value + self.half


def _mask(prose: str) -> str:
    """Blank everything that carries digits without stating a quantity. Length is preserved."""
    chars = list(prose)

    def blank(start: int, end: int) -> None:
        for index in range(start, end):
            if chars[index] != "\n":
                chars[index] = " "

    for link in markdown.links(prose):
        close = prose.rfind("](", link.start, link.end)
        if close != -1:
            blank(close + 1, link.end)                     # keep the visible text, drop the target
    for pattern in (_URL, _HTML, _HEX_RUN, _ISO_DATE, _SLASH_DATE, _CLOCK, _SEXAGESIMAL):
        for match in pattern.finditer("".join(chars)):
            blank(match.start(), match.end())
    return "".join(chars)


def _exponent(match: re.Match) -> int | None:
    raw = match.group("e1") or match.group("e3")
    if raw is None and match.group("e2"):
        raw = match.group("e2").translate(SUPERSCRIPT)
    if raw is None:
        return None
    return int(raw.replace(MINUS, "-"))


def _table_scales(prose: str) -> dict[int, dict[int, tuple[int, bool]]]:
    """For Markdown tables: line number -> column -> (power of ten, percent) declared in the header."""
    lines = prose.splitlines()
    scales: dict[int, dict[int, tuple[int, bool]]] = {}
    index = 0
    while index < len(lines) - 1:
        if _TABLE_ROW.match(lines[index]) and _TABLE_RULE.match(lines[index + 1]):
            header = [cell for cell in lines[index].strip().strip("|").split("|")]
            columns: dict[int, tuple[int, bool]] = {}
            for column, cell in enumerate(header):
                power = 0
                found = _TEN_TO.search(cell)
                if found:
                    raw = found.group(1) or found.group(2).translate(SUPERSCRIPT)
                    power = int(raw.replace(MINUS, "-"))
                percent = "%" in cell or "percent" in cell.lower()
                if power or percent:
                    columns[column] = (power, percent)
            row = index + 2
            while row < len(lines) and _TABLE_ROW.match(lines[row]):
                if columns:
                    scales[row + 1] = columns
                row += 1
            index = row
        else:
            index += 1
    return scales


def _column_of(line_text: str, offset_in_line: int) -> int:
    """Column index of a table cell, counting unescaped pipes before the offset."""
    head = line_text[:offset_in_line]
    pipes = len(re.findall(r"(?<!\\)\|", head))
    return pipes - 1 if line_text.lstrip().startswith("|") else pipes


def extract(text: str) -> list[Number]:
    """Every number stated in a Markdown document, in order of appearance."""
    prose, _ = markdown.split_code(text)
    masked = _mask(prose)
    lines = markdown.Lines(masked)
    line_starts = [0]
    for found in re.finditer("\n", masked):
        line_starts.append(found.end())
    table_scales = _table_scales(masked)
    numbers: list[Number] = []
    previous_end = -1            # end offset of the previous numeric token, kept or not
    previous_kept = True
    previous_identifier = False
    position = 0
    while True:
        match = _NUMBER.search(masked, position)
        if not match:
            break
        start, end = match.start(), match.end()
        position = end
        sign = match.group("sign") or ""
        digits_start = start + len(sign)
        whole = match.group("int")
        fraction = match.group("frac") or ""
        exponent = _exponent(match)

        # A bare "10" followed by an exponent is a power of ten: 10⁻⁸, 10^{-12}.
        bare_power = None
        if exponent is None and whole == "10" and not fraction:
            power = _POWER_OF_TEN.match(masked, end)
            if power:
                raw = power.group(1) or power.group(2).translate(SUPERSCRIPT)
                bare_power = int(raw.replace(MINUS, "-"))
                end = position = power.end()
                if masked[max(0, start - 8):start].rstrip().endswith(("×", "\\times", "\\cdot", "·", "*")):
                    previous_end, previous_kept, previous_identifier = end, False, False
                    continue                                # the multiplier of a bracket, not a number of its own

        before = masked[digits_start - 1] if digits_start > 0 else " "
        before_sign = masked[start - 1] if start > 0 else " "
        after = masked[end] if end < len(masked) else " "
        after2 = masked[end + 1] if end + 1 < len(masked) else " "
        left = masked[max(0, start - 40):start]

        def skip(identifier: bool = False) -> None:
            nonlocal previous_end, previous_kept, previous_identifier
            previous_end, previous_kept, previous_identifier = end, False, identifier

        glued_to_previous = previous_end != -1 and start <= previous_end + 1 and (
            start == previous_end or masked[previous_end] in GLUE)

        # ---- left boundary: is this token part of a longer word or code?
        signed = False
        if sign:
            if (before_sign.isalnum() and not _after_latex_command(masked, start)) or before_sign in ")]}":
                if before_sign.isdigit() or before_sign in ")]}":
                    if sign == "+" or (glued_to_previous and (previous_identifier or not previous_kept)):
                        skip(identifier=True)               # J0720+6439, 4U 1538−522, GX 1+4
                        continue
                    sign = ""                               # a range or a subtraction: 0.39-0.99
                    start = digits_start
                else:
                    skip(identifier=True)                   # KOI-7865, CD−27, Kepler-47
                    continue
            else:
                signed = True
        else:
            if (before.isalnum() or before == "_") and not _after_latex_command(masked, digits_start):
                skip(identifier=True)                       # J0720, v1, K2C9
                continue
            if before in ".," and digits_start >= 2 and masked[digits_start - 2].isdigit():
                skip(identifier=not previous_kept)          # 7.6.1 or 1,2 without spaces
                continue
            if before in "^_" or (before == "{" and digits_start >= 2 and masked[digits_start - 2] in "^_"):
                skip(identifier=True)                       # LaTeX sub- and superscripts
                continue
            if before == "/" and digits_start >= 2 and masked[digits_start - 2].isalpha():
                skip(identifier=True)                       # S/2023, C/2020
                continue
            if before in "/:" and glued_to_previous and not previous_kept:
                skip(identifier=True)
                continue
            if before in "#§":
                skip()
                continue
            if before == "–" and glued_to_previous and previous_identifier:
                skip(identifier=True)
                continue

        # ---- right boundary
        if (after.isascii() and (after.isalnum() or after == "_")) or (after == "." and after2.isdigit()) \
                or (after == ":" and after2.isdigit()) or (after == "/" and after2.isalpha()) \
                or (after == "+" and after2.isdigit()):
            skip(identifier=True)                           # 2MASS, 4U, 3D, 240P/NEAT, 1+4
            continue

        digits = whole.replace(",", "") + fraction
        try:
            mantissa = float(digits)
        except ValueError:
            skip()
            continue
        decimals = len(fraction) - 1 if fraction else 0
        scale = 10.0 ** exponent if exponent is not None else 1.0
        value = mantissa * scale
        ulp = (10.0 ** -decimals) * scale
        stripped = (whole.replace(",", "") + fraction[1:]).lstrip("0")
        sig = max(1, len(stripped))
        integer = not fraction and exponent is None and bare_power is None
        if bare_power is not None:                          # 10^n states one significant digit
            value = ulp = 10.0 ** bare_power
            sig = 1
        grouped = "," in whole
        approximate = bool(_APPROX.search(left))
        readings: list[tuple[float, float, str]] = []

        kind = "value"
        if _REFERENCE.search(left) and exponent is None and bare_power is None:
            kind = "reference"                              # Table 3, Section 5.2, version 2.1
        elif integer:
            bare = whole.replace(",", "")
            trailing = len(bare) - len(bare.rstrip("0")) if bare.strip("0") else 0
            if _MONTH_BEFORE.search(left) or _MONTH_AFTER.match(masked[end:end + 14]) \
                    or (len(bare) <= 2 and _DAY_RANGE_AFTER.match(masked[end:end + 40])):
                kind = "date"
            elif len(bare) == 4 and not grouped and 1800 <= int(bare) <= 2199:
                kind = "year"
            elif (len(bare) >= 9 and not grouped) or (len(bare) > 1 and bare.startswith("0")):
                kind = "identifier"
            elif not grouped and not signed and _DESIGNATION.search(left):
                kind = "identifier"
            elif glued_to_previous and previous_identifier:
                kind = "identifier"
            elif previous_identifier and not previous_kept and start == previous_end + 1 \
                    and masked[previous_end] == " " and not grouped:
                kind = "identifier"                         # CD−27 2812, BD+19 1234: the second half of a name
            elif before == "(" and after == ")" and len(bare) >= 3 and not grouped \
                    and re.match(r"\)\s(?:[A-Z][a-z]+|\d{4}\s[A-Z]{1,2}\d*)", masked[end:end + 12]):
                kind = "identifier"                         # (41013) 1999 UH22, (664) Judith
            if kind == "value" and trailing:
                # 205,000 may be exact or rounded to the thousand. Count only the digits that are certain.
                sig = max(1, len(bare.rstrip("0")))
                if approximate or not grouped or trailing >= 3:
                    readings.append((value, 10.0 ** trailing, "rounded to its trailing zeros"))

        if kind == "value" and value == 0 and not fraction:
            sig = 1

        if sign in ("-", MINUS):
            value = -value
        percent = masked[end:end + 2].lstrip().startswith("%") or masked[end:end + 9].lstrip().lower().startswith("percent")
        if percent:
            readings.append((value / 100.0, ulp / 100.0, "percent as a fraction"))
        multiplier = _MULTIPLIER.match(masked[end:end + 10])
        if multiplier and kind == "value":
            power = _MULTIPLIERS[multiplier.group(1).lower()]
            readings.append((value * 10.0 ** power, ulp * 10.0 ** power, f"{multiplier.group(1).lower()} applied"))

        line = lines.at(start)
        if line in table_scales and kind == "value":
            line_start = line_starts[line - 1]
            line_end = masked.find("\n", line_start)
            row_text = masked[line_start:line_end if line_end != -1 else len(masked)]
            declared = table_scales[line].get(_column_of(row_text, start - line_start))
            if declared:
                power, header_percent = declared
                if power and exponent is None:
                    readings.append((value * 10.0 ** power, ulp * 10.0 ** power, f"column header ×10^{power}"))
                if header_percent and not percent:
                    readings.append((value / 100.0, ulp / 100.0, "column header percent"))

        context = " ".join((masked[max(0, start - 36):start] + "⟦" + masked[start:end] + "⟧"
                            + masked[end:end + 28]).split())
        numbers.append(Number(text=masked[start:end], value=value, ulp=ulp, sig=sig, kind=kind, line=line,
                              start=start, end=end, percent=percent, approximate=approximate, signed=signed,
                              context=context, readings=readings))
        previous_end, previous_kept, previous_identifier = end, True, kind == "identifier"

    _share_exponents(masked, numbers)
    return numbers


def _after_latex_command(masked: str, at: int) -> bool:
    r"""True when the letters before ``at`` are a LaTeX command such as \pm or \times, not a name."""
    cursor = at - 1
    while cursor >= 0 and masked[cursor].isascii() and masked[cursor].isalpha():
        cursor -= 1
    return cursor >= 0 and cursor < at - 1 and masked[cursor] == "\\"


def _share_exponents(masked: str, numbers: list[Number]) -> None:
    """"3.85 and 4.49 × 10⁻⁸" and "(−9.8 ± 2.8) × 10⁻¹²" give every listed number the exponent."""
    inherited: dict[int, int] = {}
    for index in range(len(numbers) - 2, -1, -1):
        here, following = numbers[index], numbers[index + 1]
        if here.kind != "value" or following.kind != "value":
            continue
        if _has_exponent(here.text) or not _CONNECTOR.match(masked[here.end:following.start]):
            continue
        power = _printed_power(following)
        if power is None:
            power = inherited.get(index + 1)                # 1.2, 3.4 and 5.6 × 10⁻³
        if power is None:
            continue
        inherited[index] = power
        factor = 10.0 ** power
        here.readings.append((here.value * factor, here.ulp * factor, f"exponent ×10^{power} shared with the next number"))
    for group in re.finditer(r"\(([^()\n]{1,80})\)\s?(?:×|\\times|\\cdot|·)\s?10(?:\^\{?\s?([+\-−]?\d+)\s?\}?|([⁺⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))",
                             masked):
        raw = group.group(2) or group.group(3).translate(SUPERSCRIPT)
        factor = 10.0 ** int(raw.replace(MINUS, "-"))
        for number in numbers:
            if group.start(1) <= number.start and number.end <= group.end(1) and number.kind == "value":
                number.readings.append((number.value * factor, number.ulp * factor, "exponent applied to the bracket"))


def _has_exponent(text: str) -> bool:
    return bool(re.search(r"10\s?(?:\^|[⁺⁻⁰¹²³⁴⁵⁶⁷⁸⁹])|[eE][+\-−]?\d", text))


def _printed_power(number: Number) -> int | None:
    """The power of ten a number was printed with, or None for plain notation."""
    match = re.search(r"10(?:\^\{?\s?([+\-−]?\d+)\s?\}?|([⁺⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))\s*$", number.text)
    if match:
        raw = match.group(1) or match.group(2).translate(SUPERSCRIPT)
        return int(raw.replace(MINUS, "-"))
    match = re.search(r"[eE]([+\-−]?\d+)$", number.text)
    if match:
        return int(match.group(1).replace(MINUS, "-"))
    return None


def decimals_of(text: str) -> int | None:
    """Decimal places of a machine-written numeric string such as '12.990' or '1.5e-3'."""
    match = re.fullmatch(r"\s*[+\-]?(\d+)(?:\.(\d*))?(?:[eE]([+\-]?\d+))?\s*", text)
    if not match:
        return None
    places = len(match.group(2) or "")
    if match.group(3):
        places -= int(match.group(3))
    return places


def round_matches(candidate: float, value: float, half: float) -> bool:
    """True when ``candidate`` would be printed as ``value`` at the claimed precision."""
    slack = half * 1e-9 + abs(value) * 4 * math.ulp(1.0)
    return abs(candidate - value) <= half + slack
