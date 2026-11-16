"""Extract the two dated optical water labels, versioned GEDI shots, and provenance.

Usage: python prepare.py cache --fetch
Dependencies: Python with GDAL bindings and NumPy. cache holds original Sentinel SCL
files and local SlideRule binary responses, or --fetch downloads missing inputs.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import urllib.request

import numpy as np
from osgeo import gdal, osr

ROOT = Path(__file__).resolve().parent
BOX = (15.14, -4.41, 15.32, -4.27)
API = "https://sliderule.slideruleearth.io/source/gedi02a"
FIELDS = ["granule", "shot_number", "utc", "beam", "flags", "lat", "lon",
          "elev_lowestmode_m", "elev_highestreturn_m", "sensitivity", "solar_elevation_deg",
          "reach_id", "along_km_from_downstream_end", "offset_m",
          "jan_scl", "jul_scl", "jan_water_3x3_count", "jul_water_3x3_count"]
LAYOUT = struct.Struct("<QqddffffIBBH")


def get(cache_file, request, fetch):
    if not cache_file.exists() and fetch:
        with urllib.request.urlopen(request, timeout=120) as response:
            cache_file.write_bytes(response.read())
    return cache_file.read_bytes()


def records(raw):
    start = 0
    while start < len(raw):
        version, key_size, size = struct.unpack_from(">hhi", raw, start)
        assert version == 2 and key_size > 0 and size >= 0
        start += 8
        key = raw[start:start + key_size - 1].decode("ascii")
        start += key_size
        payload = raw[start:start + size]
        assert len(payload) == size
        start += size
        yield key, payload


def snapped_coordinates(lines, x, y):
    delta = np.array([x, y]) - lines["origins"]
    frac = np.clip(np.sum(delta * lines["segments"], axis=1) / lines["sq_lengths"], 0, 1)
    offset = np.linalg.norm(delta - lines["segments"] * frac[:, None], axis=1)
    k = np.argmin(offset)
    km = (lines["start_dist"][k] + frac[k] * np.sqrt(lines["sq_lengths"][k])) / 1000
    return int(lines["identifiers"][k]), float(km), float(offset[k])


def prepare_geometries():
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)
    source.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    target = osr.SpatialReference()
    target.ImportFromEPSG(32733)
    transform = osr.CoordinateTransformation(source, target)
    pieces = {k: [] for k in ["origins", "segments", "sq_lengths", "start_dist", "identifiers"]}
    for feature in json.loads((ROOT / "centerlines.geojson").read_text())["features"]:
        points = np.array(transform.TransformPoints(feature["geometry"]["coordinates"]))[:, :2]
        segments = np.diff(points, axis=0)
        sq_lengths = np.sum(segments * segments, axis=1)
        nonzero = sq_lengths > 0
        origins = points[:-1][nonzero]
        segments = segments[nonzero]
        sq_lengths = sq_lengths[nonzero]
        start_dist = np.r_[0, np.cumsum(np.sqrt(sq_lengths))[:-1]]
        ids = np.full(len(segments), feature["properties"]["reach_id"])
        for k, v in zip(pieces, (origins, segments, sq_lengths, start_dist, ids)):
            pieces[k].append(v)
    return transform, {k: np.concatenate(segs) for k, segs in pieces.items()}


def load_optical(cache, fetch):
    scenes = json.loads((ROOT / "optical-sources.json").read_text())["scenes"]
    rasters = []
    for scene in scenes:
        path = cache / (scene["id"] + ".tif")
        raw = get(path, scene["url"], fetch)
        assert hashlib.sha256(raw).hexdigest() == scene["sha256"]
        dataset = gdal.Open(str(path))
        inverse = gdal.InvGeoTransform(dataset.GetGeoTransform())
        rasters.append((dataset.ReadAsArray(), inverse))
    return rasters


def optical_labels(rasters, x, y):
    centers, widths = [], []
    for arr, inverse in rasters:
        col, row = map(math.floor, gdal.ApplyGeoTransform(inverse, x, y))
        if not (0 < row < arr.shape[0] - 1 and 0 < col < arr.shape[1] - 1):
            centers.append(0)
            widths.append(0)
            continue
        centers.append(int(arr[row, col]))
        widths.append(int(np.sum(arr[row-1:row+2, col-1:col+2] == 6)))
    return centers + widths


def main(cache, fetch):
    gdal.UseExceptions()
    osr.UseExceptions()
    cache.mkdir(exist_ok=True, parents=True)
    transform, lines = prepare_geometries()
    rasters = load_optical(cache, fetch)
    west, south, east, north = BOX
    polygon = [{"lon": lon, "lat": lat} for lon, lat in
               [(west, south), (east, south), (east, north), (west, north), (west, south)]]
    items = json.loads((ROOT / "cmr-granules.json").read_text())["feed"]["entry"]
    manifest = []
    # Fix gzip timestamp so repeated runs give identical derived CSV bytes.
    with open(ROOT / "gedi-corridor-shots.csv.gz", "wb") as compressed, \
            gzip.GzipFile(fileobj=compressed, mode="wb", mtime=0) as binary, \
            io.TextIOWrapper(binary, encoding="utf-8", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(FIELDS)
        for item in items:
            name = item["title"] + ".h5"
            request = {"resource": name, "parms": {"asset": "gedil2a", "poly": polygon}}
            req = urllib.request.Request(API, data=json.dumps(request).encode(),
                                         headers={"Content-Type": "application/json"})
            raw = get(cache / (name + ".sr"), req, fetch)
            n_all = 0
            missing = set()
            rows = []
            for key, payload in records(raw):
                if key == "exceptrec":
                    msg = payload[8:].split(b"\0")[0].decode(errors="replace")
                    if "Failure on resource" in msg:
                        missing.add(msg.split(" beam ")[-1].split(":")[0])
                if key != "gedi02arec":
                    continue
                for shot, ns, lat, lon, z, h, solar, sens, orbit, beam, flags, track in LAYOUT.iter_unpack(payload):
                    n_all += 1
                    x, y, _ = transform.TransformPoint(lon, lat)
                    reach, km, offset = snapped_coordinates(lines, x, y)
                    if offset > 1500:
                        continue
                    seconds, nsec = divmod(ns, 10**9)
                    utc = datetime.datetime.fromtimestamp(seconds, datetime.timezone.utc)
                    stamp = f"{utc:%Y-%m-%dT%H:%M:%S}.{nsec:09d}Z"
                    rows.append((name, shot, stamp, beam, flags, lat, lon, z, h, sens,
                                 solar, reach, round(km, 5), round(offset, 2), *optical_labels(rasters, x, y)))
            writer.writerows(sorted(rows, key=lambda r: r[1]))
            manifest.append((name, item["id"], item["time_start"], item["time_end"],
                             n_all, len(rows), ",".join(sorted(missing)), len(raw),
                             hashlib.sha256(raw).hexdigest()))
    with (ROOT / "granule-manifest.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["granule", "cmr_concept_id", "granule_start", "granule_end",
                         "shots_returned", "within_1500m", "reported_beam_read_failures",
                         "response_bytes", "local_response_sha256"])
        writer.writerows(manifest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    main(args.cache, args.fetch)
