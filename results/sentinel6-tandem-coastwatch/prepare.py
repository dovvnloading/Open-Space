"""Build a compact, replayable regional extract from the frozen source manifest."""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

from netCDF4 import Dataset
import numpy as np
import shapefile

ROOT = Path(__file__).resolve().parent
EPOCH = datetime(1985, 1, 1, tzinfo=timezone.utc)
EARTH_KM = 6371.0
STATIONS = (
    "46011", "46013", "46014", "46015", "46022", "46025", "46026", "46027",
    "46028", "46029", "46041", "46042", "46050", "46053", "46054", "46069",
)
FIELDS = {
    "time": "seconds since 1985-01-01 00:00:00 UTC",
    "lat": "degrees_north",
    "lon": "degrees_east",
    "sla": "m",
    "swh": "m",
}


def read_inputs(root):
    manifest = json.loads((ROOT / "sources.json").read_bytes())
    files = {}
    for item in manifest["indexes"] + manifest["files"]:
        if item["status"] != 200:
            continue
        source = root / item["file"]
        content = source.read_bytes()
        if sha256(content).hexdigest() != item["sha256"]:
            raise ValueError(f"Digest mismatch for {source}")
        files[item["file"]] = source
    return files


def read_satellite(path):
    with Dataset(path) as source:
        if source.featureType != "trajectory" or source.filename != path.name:
            raise ValueError(f"Unexpected daily product header in {path}")
        cols = []
        for name, units in FIELDS.items():
            variable = source.variables[name]
            if variable.units != units or variable.dimensions != ("time",):
                raise ValueError(f"Unexpected {name} units or shape in {path}")
            cols.append(np.ma.asarray(variable[:], dtype=float).filled(np.nan))
        time, lat, lon, _, _ = cols
        lon = (lon + 180.0) % 360.0 - 180.0
        cols[2] = lon
        region = np.isfinite(time) & (lat >= 32) & (lat <= 49) & (lon >= -132) & (lon <= -116)
        series = np.column_stack(cols)[region]
        metadata = {
            "native_rows": len(time),
            "regional_rows": len(series),
            "mission": source.mission_name,
            "phase": source.mission_phase,
            "history": source.history,
            "sla_quality_flag_names": source.variables["sla"].quality_flag,
            "provided_fields": sorted(source.variables),
            "sla_scale_factor_to_m": float(source.variables["sla"].scale_factor),
            "swh_scale_factor_to_m": float(source.variables["swh"].scale_factor),
        }
    return series, metadata


def coastal_points(path):
    points = []
    with ZipFile(path) as image:
        filename = next(name for name in image.namelist() if name.endswith(".shp"))
        vector = shapefile.Reader(shp=BytesIO(image.read(filename)))
        for shape in vector.iterShapes():
            if shape.bbox[0] > -111 or shape.bbox[2] < -137 or shape.bbox[1] > 54 or shape.bbox[3] < 27:
                continue
            for lower, upper in zip(shape.parts, [*shape.parts[1:], len(shape.points)], strict=True):
                line = np.array(shape.points[lower:upper], dtype=float)
                for start, stop in zip(line[:-1], line[1:], strict=True):
                    if (max(start[0], stop[0]) < -137 or min(start[0], stop[0]) > -111
                            or max(start[1], stop[1]) < 27 or min(start[1], stop[1]) > 54):
                        continue
                    # Each degree is at most 112 km on the spherical metric.
                    steps = max(1, int(np.ceil(np.linalg.norm(stop - start) * 112)))
                    points.extend(np.linspace(start[::-1], stop[::-1], steps + 1))
    if not points:
        raise ValueError("The frozen shoreline provides no regional segments")
    return np.unique(points, axis=0)


def read_buoys(files):
    locations = ElementTree.parse(files["station-locations.xml"]).getroot()
    coordinates = {
        node.attrib["id"]: (float(node.attrib["lat"]), float(node.attrib["lon"]))
        for node in locations.iter("station") if "lat" in node.attrib and "lon" in node.attrib
    }
    station_latlon = np.array([coordinates.get(station, (np.nan, np.nan)) for station in STATIONS])
    first = (datetime(2026, 9, 8, tzinfo=timezone.utc) - EPOCH).total_seconds()
    last = (datetime(2026, 11, 17, tzinfo=timezone.utc) - EPOCH).total_seconds()
    met = {}
    metadata = {}
    for index, station in enumerate(STATIONS):
        filename = f"buoys/{station}.txt"
        metadata[station] = {"coordinates": list(coordinates[station]) if station in coordinates else None,
                             "retained_wvht_rows": 0}
        if filename not in files:
            continue
        header, units, *rows = files[filename].read_text("utf-8").splitlines()
        names = header.lstrip("#").split()
        fields = [names.index(name) for name in ("YY", "MM", "DD", "hh", "mm")]
        wave = names.index("WVHT")
        if units.lstrip("#").split()[wave] != "m":
            raise ValueError(f"Unexpected WVHT unit at {station}")
        for row in rows:
            fields_in_row = row.split()
            if not fields_in_row or fields_in_row[wave] == "MM":
                continue
            instant = datetime(*[int(fields_in_row[k]) for k in fields], tzinfo=timezone.utc)
            seconds = (instant - EPOCH).total_seconds()
            value = float(fields_in_row[wave])
            if first <= seconds < last and np.isfinite(value) and value >= 0:
                met[index, seconds] = value
        metadata[station]["retained_wvht_rows"] = sum(row[0] == index for row in met)
    buoy = np.array([(index, seconds, height) for (index, seconds), height in sorted(met.items())])
    return station_latlon, buoy.reshape(-1, 3), metadata


def prepare(data_root, output, metadata_path):
    files = read_inputs(data_root)
    satellite = {"6a": [], "6b": []}
    metadata = {}
    for product, path in files.items():
        if product.startswith("satellite/"):
            data, native = read_satellite(path)
            satellite[path.name[:2]].append(data)
            metadata[path.name] = native
    extracted = {}
    for mission, daily in satellite.items():
        series = np.concatenate(daily)
        order = np.argsort(series[:, 0], kind="stable")
        extracted[mission] = series[order]
    station_latlon, buoy, station_metadata = read_buoys(files)
    coast = coastal_points(files["ne_10m_coastline.zip"])
    np.savez_compressed(output, **extracted, coast=coast, station_ids=np.array(STATIONS),
                        station_latlon=station_latlon, buoy=buoy)
    receipt = {
        "satellite_columns": list(FIELDS),
        "satellite": metadata,
        "stations": station_metadata,
        "coastline_points": len(coast),
        "data_sha256": sha256(output.read_bytes()).hexdigest(),
    }
    metadata_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.data_root, args.output, args.metadata)
