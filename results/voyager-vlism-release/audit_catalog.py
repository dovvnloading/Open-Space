"""Reproduce Voyager catalog timestamps, coverage counts and density arithmetic."""

import argparse
import csv
from datetime import datetime
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
PROFILES = {
    "voyager_bundle": ("bundle_voyager-pws-vlism-density.lblx", "ds-view-bundle.html", "api-bundle.json"),
    "voyager_collection": ("collection.lblx", "ds-view-collection.html", "api-collection.json"),
    "trajectory_control": ("xml-control-trajectory.xml", "html-control-trajectory.html", "api-control-trajectory.json"),
    "psyche_control": ("xml-control-psyche.xml", "html-control-psyche.html", "api-control-psyche.json"),
    "maven_control": ("xml-control-sep-226.xml", "html-control-sep.html", None),
}


class ProfileCells(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_cell = False
        self.text = []
        self.cells = []

    def handle_starttag(self, tag, attrs):
        if tag == "td":
            self.in_cell = True
            self.text = []

    def handle_endtag(self, tag):
        if tag == "td" and self.in_cell:
            self.cells.append(" ".join(" ".join(self.text).split()))
            self.in_cell = False

    def handle_data(self, data):
        if self.in_cell:
            self.text.append(data)


def read_xml(content, tag):
    return next(node.text.strip() for node in ET.fromstring(content).iter()
                if node.tag.endswith("}" + tag))


def utc(value):
    time = datetime.fromisoformat(value.strip())
    if time.utcoffset() is None or time.utcoffset().total_seconds() != 0:
        raise ValueError(f"Expected timezone-aware UTC: {value}")
    return time


def format_utc(time):
    return time.isoformat(timespec="milliseconds").replace("+00:00", "Z")


def label_id(label):
    return read_xml(label, "logical_identifier") + "::" + read_xml(label, "version_id")


def label_bounds(label):
    return [utc(read_xml(label, tag)) for tag in ("start_date_time", "stop_date_time")]


def api_bounds(api):
    return [utc(api[tag]) for tag in ("start_date_time", "stop_date_time")]


def body_entries(html):
    parser = ProfileCells()
    parser.feed(html.decode("utf-8", errors="replace"))
    result = {}
    for label in ("IDENTIFIER", "START DATE TIME", "STOP DATE TIME"):
        result[label] = parser.cells[parser.cells.index(label) + 1]
    return result


def archive_members():
    expected = json.loads((ROOT / "catalog-sources.json").read_text(encoding="utf-8"))["members"]
    with ZipFile(ROOT / "catalog-replay.zip") as bundle:
        results = {}
        for item in expected:
            body = bundle.read(item["member"])
            if len(body) != item["bytes"] or hashlib.sha256(body).hexdigest() != item["archive_sha256"]:
                raise ValueError(f"Recorded archive member changed: {item['member']}")
            results[item["member"]] = body
    return results


def product_rows(files, spacecraft):
    filename = f"vg{spacecraft}.csv"
    label = files[f"vg{spacecraft}.lblx"]
    data = files[filename]
    if len(data) != int(read_xml(label, "file_size")) or hashlib.md5(data).hexdigest() != read_xml(
        label, "md5_checksum"
    ):
        raise ValueError(f"Source PDS label checksum or size mismatch for {filename}")
    parsed = [
        {key.strip().split()[0]: value.strip() for key, value in raw.items()}
        for raw in csv.DictReader(io.StringIO(data.decode("ascii")))
    ]
    if len(parsed) != int(read_xml(label, "records")):
        raise ValueError(f"Source PDS row count mismatch for {filename}")
    return parsed


def read_profile(files, paths):
    label, page, api = paths
    ld = files[label]
    ids = [label_id(ld)]
    native = label_bounds(ld)
    rendered = body_entries(files[page])
    ids.append(rendered["IDENTIFIER"])
    displayed = [utc(rendered[tag]) for tag in ("START DATE TIME", "STOP DATE TIME")]
    registry = None
    if api is not None:
        api_data = json.loads(files[api])
        ids.append(api_data["id"])
        registry = api_bounds(api_data)
    matched = len(set(ids)) == 1
    folded = [time.replace(hour=time.hour % 12 or 12) for time in native]
    return {
        "ids": ids, "matched_product_version": matched,
        "source_label": [format_utc(t) for t in native],
        "public_display": [format_utc(t) for t in displayed],
        "registry_api": [format_utc(t) for t in registry] if registry else None,
        "registry_matches_label": native == registry if registry else None,
        "display_minus_label_hours": [(shown - reference).total_seconds() / 3600
                                       for shown, reference in zip(displayed, native)] if matched else None,
        "display_matches_12h_without_meridiem": displayed == folded if matched else None,
    }


def missing_rows(times, bounds):
    early = [format_utc(time) for time in times if time < bounds[0]]
    late = [format_utc(time) for time in times if time > bounds[1]]
    return {"rows_before": early, "rows_after": late, "excluded": len(early) + len(late)}


def density_arithmetic(rows):
    pairs = [(float(row["N_e"]), float(row["Fpe"])) for row in rows]
    if any(not (density > 0 and frequency > 0) for density, frequency in pairs):
        raise ValueError("Unusable archive density or plasma frequency")
    theoretical = [(frequency / 8980) ** 2 for density, frequency in pairs]
    differences = [abs(density - forecast) for (density, _), forecast in zip(pairs, theoretical)]
    by_source = {}
    for i, row in enumerate(rows):
        entry = by_source.setdefault(row["Source"].lower(), {
            "rows": 0, "max_absolute_difference_cm3": 0., "max_relative_difference": 0.,
            "invalid_density_bounds": 0,
        })
        measured, _ = pairs[i]
        entry["rows"] += 1
        entry["max_absolute_difference_cm3"] = max(entry["max_absolute_difference_cm3"], differences[i])
        entry["max_relative_difference"] = max(entry["max_relative_difference"], differences[i] / measured)
        low, high = float(row["N_e_min"]), float(row["N_e_max"])
        entry["invalid_density_bounds"] += int(not 0 < low <= measured <= high)
    return by_source


def run():
    files = archive_members()
    comparisons = {name: read_profile(files, paths) for name, paths in PROFILES.items()}
    vprofile, cprofile = [comparisons[k] for k in ("voyager_bundle", "voyager_collection")]
    cover = {
        "bundle_xml": vprofile["source_label"], "bundle_display": vprofile["public_display"],
        "collection_xml": cprofile["source_label"], "collection_display": cprofile["public_display"],
    }
    products = {}
    for craft in (1, 2):
        items = product_rows(files, craft)
        times = [utc(record["SCET"]) for record in items]
        label = files[f"vg{craft}.lblx"]
        api = json.loads(files[f"api-voyager{craft}.json"])
        if label_id(label) != api["id"]:
            raise ValueError("The product registry returned a different version")
        limits = [min(times), max(times)]
        lidvid = label_id(label)
        inventory = list(csv.reader(io.StringIO(files["inventory.csv"].decode("ascii"))))
        products[f"voyager{craft}"] = {
            "lidvid": lidvid, "rows": len(items), "unique_times": len(set(times)),
            "first_last_scet": [format_utc(t) for t in limits],
            "listed_in_collection_inventory": any(lidvid in line for line in inventory),
            "scet_matches_product_label": limits == label_bounds(label),
            "scet_matches_registry": limits == api_bounds(api),
            "filter_effects": {name: missing_rows(times, [utc(t) for t in limits])
                               for name, limits in cover.items()},
            "density_conversion": density_arithmetic(items),
        }
    return {
        "exact_version_profile_comparisons": comparisons,
        "registry_api_unavailable_for_matching_maven_version": True,
        "products": products,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n",
                           encoding="utf-8", newline="\n")
