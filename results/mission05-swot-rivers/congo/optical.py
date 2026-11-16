"""Measure Sentinel-2B water on identical clear pixels around the SWORD v16 centerline.

Run with Python 3, NumPy and the GDAL/OGR Python bindings. The script verifies the
original Sentinel COG SHA-256 hashes and recreates the small archived three-band crop.
"""

import hashlib
import json
import math
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
from osgeo import gdal, ogr, osr

ROOT = Path(__file__).resolve().parent


def read_source(directory, source):
    path = directory / (source["id"] + ".tif")
    with urllib.request.urlopen(source["url"], timeout=60) as response:
        data = response.read()
    assert hashlib.sha256(data).hexdigest() == source["sha256"]
    path.write_bytes(data)
    return gdal.Open(str(path))


def measure(raster):
    jan, jul, mask = raster.ReadAsArray()
    roi = mask.astype(bool)
    common = roi & np.isin(jan, [4, 5, 6]) & np.isin(jul, [4, 5, 6])
    return {
        "corridor_pixels": int(roi.sum()),
        "common_clear_pixels": int(common.sum()),
        "high_water_pixels": int((common & (jan == 6)).sum()),
        "low_water_pixels": int((common & (jul == 6)).sum()),
        "pixel_area_m2": 400,
    }


def main():
    gdal.UseExceptions()
    ogr.UseExceptions()
    config = json.loads((ROOT / "independent/optical-source.json").read_text())
    geometry = ogr.Open(str(ROOT / "sword-13217000061.geojson"))
    layer = geometry.GetLayer()
    feature = layer.GetNextFeature()
    line = feature.GetGeometryRef().Clone()
    target = osr.SpatialReference()
    target.ImportFromEPSG(32733)
    target.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    source = layer.GetSpatialRef()
    source.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    line.Transform(osr.CoordinateTransformation(source, target))
    corridor = line.Buffer(1000)
    xmin, xmax, ymin, ymax = corridor.GetEnvelope()
    with tempfile.TemporaryDirectory() as temp:
        images = [read_source(Path(temp), item) for item in config["scenes"]]
        gt = images[0].GetGeoTransform()
        assert images[1].GetGeoTransform() == gt
        c = math.floor((xmin - gt[0]) / 20)
        r = math.floor((gt[3] - ymax) / 20)
        w = math.ceil((xmax - gt[0]) / 20) - c
        h = math.ceil((gt[3] - ymin) / 20) - r
        out = gdal.GetDriverByName("GTiff").Create(
            str(ROOT / "independent/sentinel2-common-grid.tif"),
            w,
            h,
            3,
            gdal.GDT_Byte,
            options=["COMPRESS=DEFLATE"],
        )
        out.SetGeoTransform((gt[0] + c * 20, 20, 0, gt[3] - r * 20, 0, -20))
        out.SetProjection(target.ExportToWkt())
        for i, img in enumerate(images, 1):
            out.GetRasterBand(i).WriteArray(img.ReadAsArray(c, r, w, h))
            out.GetRasterBand(i).SetDescription(config["scenes"][i - 1]["id"])
        mem = ogr.GetDriverByName("Memory").CreateDataSource("")
        sheet = mem.CreateLayer("mask", target, ogr.wkbPolygon)
        cell = ogr.Feature(sheet.GetLayerDefn())
        cell.SetGeometry(corridor)
        sheet.CreateFeature(cell)
        gdal.RasterizeLayer(out, [3], sheet, burn_values=[1])
        out.GetRasterBand(3).SetDescription(
            "SWORD v16 reach 13217000061 buffered by 1000 m"
        )
        result = measure(out)
        out = None
    (ROOT / "independent/optical-area.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(result)


if __name__ == "__main__":
    main()
