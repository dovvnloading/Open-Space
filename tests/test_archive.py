import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from open_space import archive

DATA = Path(__file__).resolve().parents[1] / "data"
PILOT_NAME = "sdss-6497-56329-0139"


class ArchiveTests(unittest.TestCase):
    def test_pinned_observation_and_real_known_lines(self):
        manifest = json.loads((DATA / f"{PILOT_NAME}.json").read_text())
        self.assertEqual(manifest["source_url"], archive.archive_source(*archive.PILOT)[1])
        self.assertEqual(manifest["source_sha256"], archive.PILOT_SHA256)
        raw = (DATA / manifest["csv_file"]).read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), manifest["csv_sha256"])
        rows = list(csv.DictReader(io.StringIO(raw.decode("ascii"))))
        self.assertEqual(len(rows), manifest["pixels"])
        self.assertTrue(all(row["mask"] == row["and_mask"] for row in rows))
        usable = [r for r in rows if int(r["mask"]) == 0 and float(r["ivar"]) > 0]
        self.assertEqual(len(usable), manifest["usable_pixels"])
        # Positive controls use published rest wavelengths and survey redshift.
        for rest in (4862.683, 4960.295, 5008.239, 6564.614):
            predicted = rest * (1 + manifest["redshift"])
            window = [r for r in usable if abs(float(r["wavelength"]) - predicted) < 4]
            peak = max(window, key=lambda r: float(r["flux"]))
            self.assertLess(abs(float(peak["wavelength"]) - predicted), 1)
            self.assertGreater(int(peak["or_mask"]), 0)

    def test_source_checksum_and_output_conflicts(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "sample"
            archive.write_immutable(output, b"one")
            archive.write_immutable(output, b"one")
            with self.assertRaisesRegex(ValueError, "Conflicting output"):
                archive.write_immutable(output, b"two")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                archive.read_source("", output, "0" * 64)
            self.assertEqual(output.read_bytes(), b"one")

    def test_invalid_identifiers(self):
        for fields in [(0, 56329, 139), (6497, 999, 139), (6497, 56329, 1001)]:
            with self.assertRaises(ValueError):
                archive.archive_source(*fields)

    def test_fits_conversion_preserves_masks_and_uses_galaxy_fit(self):
        try:
            from astropy.io import fits
        except ModuleNotFoundError:
            self.skipTest("Astropy is an optional archive dependency")
        primary = fits.PrimaryHDU()
        for key, value in zip(("PLATEID", "MJD", "FIBERID"), archive.PILOT):
            primary.header[key] = value
        primary.header["VACUUM"] = True
        primary.header["RUN2D"] = archive.RUN2D
        def column(name, fmt, values):
            return fits.Column(name=name, format=fmt, array=values)
        coadd = fits.BinTableHDU.from_columns([
            column("loglam", "E", [3.5, 3.6, 3.7]),
            column("flux", "E", [5.0, float("nan"), 6.0]),
            column("ivar", "E", [2.0, 1.0, 3.0]),
            column("and_mask", "J", [0, 0, 1]),
            column("or_mask", "J", [1, 2, 3]),
        ], name="COADD")
        metadata = [
            ("PLATE", "J", 6497), ("MJD", "J", 56329), ("FIBERID", "J", 139),
            ("CLASS", "8A", "GALAXY"), ("CLASS_NOQSO", "8A", "GALAXY"),
            ("SUBCLASS_NOQSO", "10A", "STARBURST"), ("Z_NOQSO", "D", 0.1),
            ("Z_ERR_NOQSO", "D", 0.001), ("ZWARNING_NOQSO", "J", 0),
            ("PLUG_RA", "D", 1.0), ("PLUG_DEC", "D", -1.0),
            ("PLATEQUALITY", "6A", "good"), ("PROGRAMNAME", "6A", "boss"),
            ("SPECPRIMARY", "J", 1),
        ]
        spall = fits.BinTableHDU.from_columns(
            [column(name, fmt, [value]) for name, fmt, value in metadata], name="SPALL"
        )
        blob = io.BytesIO()
        fits.HDUList([primary, coadd, spall]).writeto(blob)
        raw = blob.getvalue()
        content, meta = archive.csv_and_metadata(raw, archive.PILOT)
        rows = list(csv.DictReader(io.StringIO(content.decode())))
        self.assertAlmostEqual(float(rows[0]["wavelength"]), 10**3.5)
        self.assertEqual([float(r["ivar"]) for r in rows], [2.0, 0.0, 3.0])
        self.assertEqual([r["mask"] for r in rows], ["0", "0", "1"])
        self.assertEqual([r["or_mask"] for r in rows], ["1", "2", "3"])
        self.assertEqual(meta["usable_pixels"], 1)
        self.assertEqual(meta["redshift_column"], "SPALL.Z_NOQSO")
        with self.assertRaisesRegex(ValueError, "match requested observation"):
            archive.csv_and_metadata(raw, (6497, 56329, 140))
