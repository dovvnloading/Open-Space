"""Read archived MaRS fixed-width L2 Doppler tables using their PDS3 labels."""

from datetime import datetime
import math
from pathlib import Path
import re


def _scalar(label, name):
    match = re.search(r"^\s*" + re.escape(name) + r"\s*=\s*([^\r\n]+)", label, re.M)
    if not match:
        raise ValueError(f"Missing label field {name}")
    return match.group(1).strip().strip('"')


def _columns(label):
    columns = {}
    for body in re.findall(r"^\s*OBJECT\s*=\s*COLUMN\s*$(.*?)^\s*END_OBJECT\s*=\s*COLUMN\s*$",
                           label, re.M | re.S):
        name_match = re.search(r'\bNAME\s*=\s*"([^"]+)"', body, re.S)
        if not name_match:
            raise ValueError("Unnamed label column")
        name = " ".join(name_match.group(1).split()).upper()
        start, width = int(_scalar(body, "START_BYTE")), int(_scalar(body, "BYTES"))
        if name in columns or start < 1 or width < 1:
            raise ValueError(f"Invalid or duplicated column {name}")
        columns[name] = (start - 1, width)
    return columns


def _find(columns, prefix):
    matches = [value for name, value in columns.items() if name.startswith(prefix)]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {prefix} column")
    return matches[0]


def _number(value):
    value = value.strip()
    if not value:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite table value")
    # PDS labels document these missing/corrupt and closed-loop placeholders.
    if result in (-999.9, -999.999999, -99999.999999, -999999999.999999):
        return None
    return result


def parse_rows(table_path, label_path):
    """Return rows with UTC receiver time and label-defined numeric columns.

    Geometry is the label's DISTANCE field; its meaning depends on observation
    type. Frequency values are Hz. Missing PDS sentinel values become None.
    """
    table_path, label_path = Path(table_path), Path(label_path)
    label = label_path.read_text(encoding="ascii")
    if _scalar(label, "^DOPPLER_TABLE").casefold() != table_path.name.casefold():
        raise ValueError("Label table pointer does not match supplied product")
    if _scalar(label, "RECORD_TYPE") != "FIXED_LENGTH":
        raise ValueError("Expected fixed-length records")
    row_bytes = int(_scalar(label, "ROW_BYTES"))
    if row_bytes != int(_scalar(label, "RECORD_BYTES")):
        raise ValueError("Record and row lengths differ")
    row_count = int(_scalar(label, "ROWS"))
    if row_count != int(_scalar(label, "FILE_RECORDS")):
        raise ValueError("Label row and file record counts differ")
    columns = _columns(label)
    fields = {
        "native_row": _find(columns, "SAMPLE NUMBER"),
        "utc": _find(columns, "UTC TIME"),
        "geometry_km": _find(columns, "DISTANCE"),
        "observed_hz": _find(columns, "OBSERVED "),
        "predicted_hz": _find(columns, "PREDICTED "),
        "correction_hz": _find(columns, "CORRECTION OF EARTH ATMOSPHERE"),
        "residual_hz": _find(columns, "RESIDUAL CALIBRATED "),
        "differential_hz": _find(columns, "DIFFERENTIAL DOPPLER"),
        "sigma_hz": _find(columns, "SIGMA OBSERVED ANTENNA FREQUENCY"),
    }
    if any(start + width > row_bytes - 2 for start, width in fields.values()):
        raise ValueError("Column extends beyond record payload")
    if table_path.stat().st_size != row_count * row_bytes:
        raise ValueError("Table byte length does not match label")
    result = []
    with table_path.open("rb") as stream:
        for index in range(row_count):
            raw = stream.read(row_bytes)
            if len(raw) != row_bytes or raw[-2:] != b"\r\n":
                raise ValueError(f"Invalid record length or CRLF at row {index + 1}")
            slices = {key: raw[start:start + width].decode("ascii").strip()
                      for key, (start, width) in fields.items()}
            native_row = int(slices.pop("native_row"))
            if native_row != index + 1:
                raise ValueError(f"Unexpected native row number {native_row}")
            stamp = slices.pop("utc")
            datetime.fromisoformat(stamp)
            result.append({"utc": stamp + "Z", "native_row": native_row,
                           **{name: _number(value) for name, value in slices.items()}})
        if stream.read(1):
            raise ValueError("Unexpected trailing table bytes")
    return result
