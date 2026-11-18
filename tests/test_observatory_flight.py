"""Source scope, immutable plates and physical-coordinate checks for research flights."""

import hashlib
import json
from pathlib import Path
import struct
from tempfile import TemporaryDirectory
import unittest
import zlib

from observatory.flight_plates import VECTOR_PLOTS, render_plate
from observatory.research_destinations import build_destinations, pointing, publish_image, report_excerpt, scope
from open_space.verify.tree import RepoTree


def fits(values, width, height):
    cards = {"SIMPLE": "T", "BITPIX": "-32", "NAXIS": "2", "NAXIS1": width, "NAXIS2": height,
             "DATE-OBS": "'2025-05-01'"}
    header = "".join(f"{key:<8}= {value!s:>20}".ljust(80) for key, value in cards.items()) + "END".ljust(80)
    return header.ljust(2880).encode("ascii") + struct.pack(f">{len(values)}f", *values)


class ObservatoryFlightTests(unittest.TestCase):
    repo = "https://github.com/owner/research"

    def test_registered_pixels_use_native_values_and_preserve_axis_order(self):
        raw = fits([0, 1, 2, 3, 4, float("nan")], 3, 2)
        png, date = render_plate(raw)
        self.assertEqual(date, "2025-05-01")
        # The stored axis-2 upper row goes first. 1st and 99.5th
        # lower-order percentile statistics are 0 and 3 in this small sample.
        index = png.index(b"IDAT")
        size = int.from_bytes(png[index - 4:index], "big")
        self.assertEqual(zlib.decompress(png[index + 4:index + 4 + size]),
                         bytes([0, 255, 255, 0, 0, 0, 85, 170]))
        self.assertRaises(ValueError, render_plate, fits([1, 1], 2, 1))
        self.assertRaises(ValueError, render_plate, raw.replace(b"-32", b"-64", 1))

    def test_excerpts_keep_complete_opening_prose_and_accurate_lines(self):
        excerpt = "Three observations agree with the reported trend. A variable background remains possible."
        first = f"# Status\n\n{excerpt} [Report](paper.md) · [Table](table.csv)\n\n## Reproduction\n\nRun these commands."
        self.assertEqual(report_excerpt(first, 8, ""), (excerpt, 10))
        self.assertIsNone(report_excerpt(first, 1, excerpt))
        self.assertIsNone(report_excerpt("# Calibration\n\n## Reproduction\n\nRun these commands.", 1, ""))
        self.assertIsNone(report_excerpt("# Report\n\nRun these commands to reproduce all detections from the published files.", 1, ""))

    def test_selects_only_committed_raster_and_table_from_requested_section(self):
        png, _ = render_plate(fits([0, 1, 2, 3], 2, 2))
        path = "results/field/README.md"
        document = """# Paper

## Observations

![Comparison](comparison.png)

## Observations

```md
![Code example](code.png)
```

![Cross-study](../neighbor/secret.png)
![Remote](https://example.net/site.png)
![Markup](active.svg)
![Wrong bytes](text.png)
![Visit image](field.png)

| Native epoch | Signal ± local noise (mJy/beam) |
| --- | --- |
| 2020 | −0.12 ± 0.31 |
| 2021 | 0.43 ± 0.14 |

## Another study

![Wrong field](other.png)
"""
        files = {path: document.encode(), "results/field/field.png": png, "results/field/text.png": b"not an image",
                 "results/field/code.png": png, "results/neighbor/secret.png": png,
                 "results/field/comparison.png": png, "results/field/active.svg": b"<svg/>"}
        tree = RepoTree(files, commit="1" * 40)
        row = {"path": path, "id": "results-field", "topic": "Radio", "title": "Two visits",
               "url": f"{self.repo}/blob/{tree.commit}/{path}#observations-1", "summary": ""}
        reads = []
        def read(name):
            reads.append(name)
            return tree.read(name)
        with TemporaryDirectory() as tmp:
            site = Path(tmp)
            item, = build_destinations(tree, [row], [], self.repo, site / "data/catalog.json", read)
            self.assertEqual(item["id"], item["recordId"])
            self.assertEqual(item["table"]["columns"][1], "Signal ± local noise (mJy/beam)")
            self.assertEqual(item["table"]["rows"][0], ["2020", "−0.12 ± 0.31"])
            self.assertEqual(item["table"]["totalRows"], 2)
            line = document.splitlines().index("| Native epoch | Signal ± local noise (mJy/beam) |") + 1
            self.assertTrue(item["table"]["source"]["url"].endswith(f"?plain=1#L{line}"))
            self.assertEqual(reads, ["results/field/text.png", "results/field/field.png"])
            self.assertEqual((site / item["figure"]["src"]).read_bytes(), png)
            self.assertEqual(item["figure"]["sha256"], hashlib.sha256(png).hexdigest())
            self.assertNotIn("coordinates", item)
            self.assertFalse(item["featured"])
        self.assertEqual(scope(document, "unknown"), ("", 1))

    def test_coordinate_is_from_named_row_in_reported_frame(self):
        path = "results/jwst-repeat-imaging-transients/goods-s-oct2026/README.md"
        doc = ("# Targets\n\n| Site | Position (J2000 degrees) |\n|---|---|\n"
               "| 525 | 12.3, +45.6 |\n| 185 | 53.188672, −27.812554 |")
        tree = RepoTree({path: doc.encode()}, commit="a" * 40)
        args = ("jwst-repeat-imaging-transients", doc, path, tree, tree.read, self.repo)
        coord = pointing(*args)
        self.assertEqual((coord["raDeg"], coord["decDeg"], coord["frame"]), (53.188672, -27.812554, "J2000"))
        self.assertTrue(coord["source"]["url"].endswith("#L6"))
        self.assertIsNone(pointing(args[0], doc.replace("185", "999"), *args[2:]))
        self.assertIsNone(pointing(args[0], doc.replace("J2000", "pixels"), *args[2:]))

    def test_revised_vector_source_requires_a_matching_digest(self):
        path = next(iter(VECTOR_PLOTS))
        tree = RepoTree({path: b'<svg xmlns="http://www.w3.org/2000/svg"><text>Different plot</text></svg>'},
                        commit="b" * 40)
        with TemporaryDirectory() as tmp:
            self.assertIsNone(publish_image(tree, tree.read, Path(tmp), self.repo, path, "Reviewed plot"))
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_packaged_catalog_covers_every_record_and_resolves_assets(self):
        site = Path(__file__).resolve().parents[1] / "observatory"
        catalog = json.loads((site / "data/catalog.json").read_bytes())
        destinations = catalog["destinations"]
        self.assertEqual({d["recordId"] for d in destinations}, {r["id"] for r in catalog["records"]})
        self.assertEqual(len({d["id"] for d in destinations}), len(catalog["records"]))
        self.assertTrue({s["id"] for s in catalog["spotlights"]} <= {d["id"] for d in destinations if d["featured"]})
        for d in destinations:
            with self.subTest(id=d["id"]):
                for key in ("image", "figure"):
                    if key in d:
                        asset = d[key]
                        self.assertIn(catalog["revision"], asset["url"])
                        path = site / asset["src"]
                        self.assertTrue(path.is_file(), str(path))
                        if asset.get("sha256"):
                            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), asset["sha256"])
                if table := d.get("table"):
                    self.assertTrue(all(len(row) == len(table["columns"]) for row in table["rows"]))
                    self.assertTrue(0 < len(table["rows"]) <= table["totalRows"])
                    self.assertIn(catalog["revision"], table["source"]["url"])
                if coords := d.get("coordinates"):
                    self.assertTrue(0 <= coords["raDeg"] < 360 and -90 <= coords["decDeg"] <= 90)
                    self.assertIn(catalog["revision"], coords["source"]["url"])
