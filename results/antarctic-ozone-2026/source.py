"""Read NASA daily OMPS/OMI 1-degree ASCII ozone grids without masking values.

Zeros remain zeros: their missing-value meaning must be established separately
from documentation for these particular products. Rows run south to north and
columns west to east. No analysis or interpretation is performed here.
"""
from datetime import datetime
from pathlib import Path
import re


HEADER = re.compile(
    r"^\s*Day:\s*(\d{1,3})\s+([A-Za-z]{3})\s+(\d{1,2}),\s*(\d{4})\s+"
    r"(OMPS|OMI)\s+TO3\s+STD\s+OZONE\s+GEN:\s*\d{2}:\d{3}\s+Asc\s+LECT:\s*\d{1,2}:\d{2}\s*[APap][Mm]\s*$"
)
LON_HEADER = re.compile(
    r"^\s*Longitudes:\s*360 bins centered on\s+179\.5\s+W\s+to\s+179\.5\s+E\s+\(1\.00 degree steps\)\s*$"
)
LAT_HEADER = re.compile(
    r"^\s*Latitudes\s*:\s*180 bins centered on\s+89\.5\s+S\s+to\s+89\.5\s+N\s+\(1\.00 degree steps\)\s*$"
)
MONTHS = {name: i for i, name in enumerate(
    ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 1)}


def read_grid(path, expected_date=None, expected_sensor=None):
    """Return metadata and a 180-by-360 integer list in Dobson units.

    ``expected_date`` is an optional ISO calendar date; ``expected_sensor``
    is either ``OMPS`` or ``OMI``. Reject altered grid descriptions, dates,
    latitude ordering, record counts, and malformed three-character cells.
    """
    lines = Path(path).read_text(encoding="ascii").splitlines()
    if len(lines) != 2703:
        raise ValueError(f"Expected 2703 lines; got {len(lines)}")
    match = HEADER.fullmatch(lines[0])
    if not match or not LON_HEADER.fullmatch(lines[1]) or not LAT_HEADER.fullmatch(lines[2]):
        raise ValueError("Unrecognized ozone grid header")
    doy, month, day, year, sensor = match.groups()
    if month.title() not in MONTHS:
        raise ValueError("Invalid header month")
    date = datetime(int(year), MONTHS[month.title()], int(day)).date()
    if int(doy) != date.timetuple().tm_yday:
        raise ValueError("Header day-of-year contradicts date")
    if expected_date is not None and date.isoformat() != str(expected_date):
        raise ValueError("Header date differs from expected date")
    if expected_sensor is not None and sensor != expected_sensor:
        raise ValueError("Header sensor differs from expected sensor")
    grid = []
    for band in range(180):
        row = []
        for part in range(15):
            line = lines[3 + band * 15 + part]
            count = 25 if part < 14 else 10
            width = 1 + count * 3
            if len(line) != (76 if part < 14 else 46) or line[0] != " ":
                raise ValueError(f"Malformed record at latitude band {band}, part {part}")
            body = line[1:width]
            for cell in range(count):
                text = body[cell * 3:(cell + 1) * 3]
                if not text.strip().isdigit():
                    raise ValueError(f"Invalid ozone cell at band {band}, part {part}")
                row.append(int(text))
            trailer = line[width:]
            if part == 14:
                label = re.fullmatch(r"\s+lat\s*=\s*([+-]?\d+\.\d)\s*", trailer)
                if not label or abs(float(label.group(1)) - (-89.5 + band)) > 1e-6:
                    raise ValueError(f"Unexpected latitude annotation at band {band}")
            elif trailer:
                raise ValueError("Unexpected text after ozone record")
        if len(row) != 360:
            raise ValueError("Unexpected longitude count")
        grid.append(row)
    return {"date": date.isoformat(), "sensor": sensor, "units": "DU",
            "latitude_centers": [-89.5 + i for i in range(180)],
            "longitude_centers": [-179.5 + i for i in range(360)],
            "zero_semantics": "unresolved; zeros retained"}, grid
