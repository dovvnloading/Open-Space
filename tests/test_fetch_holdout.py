import csv
import importlib.util
import itertools
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from examples import fetch_holdout as intake


def record(plate, mjd, ra, dec=10, fiber=None):
    row = {"plate": str(plate), "mjd": str(mjd), "ra": str(ra), "dec": str(dec)}
    if fiber is not None:
        row["fiber"] = str(fiber)
    return row


class HoldoutIntakeTests(unittest.TestCase):
    def test_predeclared_geometry_order_and_coverage(self):
        first = record(6497, 56329, 0)
        fields = [first, record(6498, 56500, 2), record(6499, 56501, 30),
                  record(6500, 56502, 60), record(6501, 56501, 90),
                  record(6502, 56503, 90), record(6503, 56504, 120)]
        targets = []
        for field, count in zip(fields, (8, 8, 10, 5, 8, 7, 8), strict=True):
            for i in range(count):
                row = {**field, "fiber": str(i + 1), "ra": str(float(field["ra"]) + i / 100)}
                targets.append(row)
        # The second fiber of plate 6499 is the same sky target as the first.
        target = next(r for r in targets if r["plate"] == "6499" and r["fiber"] == "2")
        target["ra"] = "30.0001"
        selected = intake.select_fields(reversed(fields), reversed(targets))
        self.assertEqual([int(r["plate"]) for r in selected], [6499]*8+[6502]*7+[6503]*8)
        self.assertEqual([int(r["fiber"]) for r in selected[:8]], [1, 3, 4, 5, 6, 7, 8, 9])
        self.assertAlmostEqual(intake.distance((359, 0), (1, 0)), 2.0)
        with self.assertRaises(ValueError):
            intake.select_fields(fields[:4], targets)

    def test_pinned_sample_and_science_exposures(self):
        with intake.TARGETS.open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 24)
        self.assertEqual(len({(r["plate"], r["mjd"], r["fiber"]) for r in rows}), 24)
        self.assertEqual(hashlib.sha256(intake.selection_csv(rows)).hexdigest(),
                         "5ec2510cdb1f69b879706995cc8115ff910396a0f5d3661853a01a2f3b3626ae")
        ids = {}
        for row in rows:
            ids.setdefault(intake.visit(row), set()).update(row["science_exposures"].split())
            self.assertRegex(row["source_sha256"], r"^[a-f0-9]{64}$")
            self.assertRegex(row["csv_sha256"], r"^[a-f0-9]{64}$")
            self.assertRegex(row["manifest_sha256"], r"^[a-f0-9]{64}$")
        self.assertEqual(len(ids), 3)
        for a, b in itertools.combinations(ids.values(), 2):
            self.assertFalse(a & b)

    @unittest.skipUnless(importlib.util.find_spec("astropy"), "Astropy required to read FITS headers")
    def test_camera_deduplication_and_invalid_header(self):
        from astropy.io import fits
        import io
        hdu = fits.PrimaryHDU()
        hdu.header["NEXP"] = 2
        hdu.header["EXPID001"] = "b1-00166857-00166855-00166856"
        hdu.header["EXPID002"] = "r1-00166857-00166855-00166856"
        raw = io.BytesIO()
        hdu.writeto(raw)
        self.assertEqual(intake.science_exposures(raw.getvalue()), ["00166857"])
        hdu.header["EXPID001"] = "unrecognized"
        raw = io.BytesIO()
        hdu.writeto(raw)
        with self.assertRaisesRegex(ValueError, "EXPID"):
            intake.science_exposures(raw.getvalue())

    def test_explicit_cache_miss_is_offline_and_preserves_data(self):
        with intake.TARGETS.open() as stream:
            first = next(csv.DictReader(stream))
        with tempfile.TemporaryDirectory() as temp, patch("open_space.archive.urlopen") as remote:
            root = Path(temp)
            with self.assertRaises(FileNotFoundError):
                intake.fetch_pinned([first], root / "out", root / "absent")
            remote.assert_not_called()
