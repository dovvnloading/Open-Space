"""Extract the published CTX time brackets and matching PDS footprints.

Usage: python extract.py CATALOG.xlsx COVERAGE.zip > endpoints.csv
Requires openpyxl and pyshp. All distances are calculated by audit.py.
"""

import csv
import json
import sys
from itertools import islice

import openpyxl
import shapefile


def endpoints(catalog):
    with open(catalog, "rb") as handle:
        workbook = openpyxl.load_workbook(handle, read_only=True, data_only=True)
        try:
            for number, row in enumerate(islice(workbook.active.values, 1, None), 2):
                for column, stage in [(13, "before"), (16, "after")]:
                    if row[column] == "CTX":
                        yield [number, row[0], row[1], row[2], stage,
                               row[column + 1].strip(), row[column + 2].date().isoformat()]
        finally:
            workbook.close()


def main(catalog, coverage):
    observations = list(endpoints(catalog))
    # Three spreadsheet strings truncate IDs. Require a unique prefix match.
    prefixes = {r[5].removeprefix("CTX ") for r in observations}
    lengths = {len(prefix) for prefix in prefixes}
    matches = {}
    with shapefile.Reader(coverage) as layer:
        for item in layer.iterRecords(fields=["ProductId", "UTCstart"]):
            name, timestamp = [field.strip() for field in item]
            for size in lengths:
                prefix = name[:size]
                if prefix not in prefixes:
                    continue
                if prefix in matches and matches[prefix][:2] != [name, timestamp]:
                    raise ValueError(f"Ambiguous prefix {prefix}")
                shape = layer.shape(item.oid)
                offsets = [*shape.parts, len(shape.points)]
                rings = [shape.points[i:j] for i, j in zip(offsets, offsets[1:])]
                matches.setdefault(prefix, [name, timestamp, []])[2].extend(rings)
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["row", "site", "lat", "lon", "stage", "catalog_id", "catalog_date",
                     "pds_product", "pds_start_utc", "footprint_rings"])
    for row in observations:
        key = row[5].removeprefix("CTX ")
        if key not in matches:
            raise ValueError(f"Missing CTX product {key}")
        name, timestamp, rings = matches[key]
        writer.writerow([*row, name, timestamp, json.dumps(rings, separators=(",", ":"))])


if __name__ == "__main__":
    main(*sys.argv[1:])
