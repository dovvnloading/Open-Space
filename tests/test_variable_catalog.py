import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from open_space import variable_catalog as cat


class VariableCatalogTest(unittest.TestCase):
    def test_validate_empty_vs_broken_tables(self):
        text = (
            b"#INFO\trequest_date=2026-01-01T00:00:00\tUTC\n#Name: A/B\n"
            b"#Title: Test\n#Column\tName\t(c)\n#Column\tPeriod\t(d)\n"
        )
        cols = "Name,Period"
        self.assertEqual(cat._parse_vizier(text, "A/B", cols)[0], [])
        full = text + b"Name\tPeriod\n \td\n---\t-----\n star\t1.50\n"
        self.assertEqual(cat._parse_vizier(full, "A/B", cols)[0], [{"Name": "star", "Period": "1.50"}])
        for broken in (
            full + b"#INFO\tError=failed",
            full.replace(b"#Name: A/B", b"#Name: A"),
            full.replace(b"#Column\tPeriod", b"#Column\twrong"),
            full.replace(b" star\t1.50", b" star\t1.50\t2.3"),
            b"<html>",
        ):
            with self.assertRaises(ValueError):
                cat._parse_vizier(broken, "A/B", cols)

    def test_query_errors_are_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(
            cat, "retrieve", side_effect=TimeoutError
        ):
            result = cat.search(180, 45, tmp)
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(result["vsx"]["status"], "error")
        self.assertEqual(result["tars"]["status"], "requires_tic_id")

    def test_exact_id_lookup_and_cache(self):
        cols = cat.RELEASES["tess_svc"][4]
        dataset = (cols + "\n42,0,-30,2sin,3,2,.02,1,.12\n101,5,-20,1sin,1,3,.2,nan,nan\n").encode()
        digest = hashlib.sha256(dataset).hexdigest()
        spec = {"tess_svc": ("https://example.test/test.csv", digest, 2, "small.csv", cols)}
        with tempfile.TemporaryDirectory() as tmp, patch.dict(cat.RELEASES, spec), patch.object(
            cat, "download", return_value=dataset
        ) as download:
            cache = Path(tmp)
            result = cat._tic_catalog("tess_svc", 42, 0.001, -30, 1.001, cache)
            self.assertEqual(result["status"], "ok")
            match = result["matches"][0]
            self.assertAlmostEqual(match["separation_arcsec"], 3.1177, places=3)
            self.assertEqual(match["periods"][0]["period_comparison"]["closest_factor"], 0.5)
            self.assertEqual(cat._tic_catalog("tess_svc", 123, 0, -30, None, cache)["matches"], [])
            download.assert_called_once()
            bad = {"tess_svc": (*spec["tess_svc"][:1], "f" * 64, *spec["tess_svc"][2:])}
            with patch.dict(cat.RELEASES, bad), self.assertRaisesRegex(ValueError, "checksum"):
                cat._index("tess_svc", cache)

    def test_known_fundamental_wider_than_discovery_range(self):
        try:
            import numpy as np
            from astropy.timeseries import LombScargle
        except ImportError:
            self.skipTest("Astropy and NumPy required")
        from open_space.stellar_harmonics import analyze_sector

        times = np.arange(0, 26, 1 / 48)
        flux = np.sin(2 * np.pi * times / 9.6) + 0.15 * np.sin(2 * np.pi * times / 4.8)
        evidence = analyze_sector(times, flux, 4.8)
        self.assertAlmostEqual(evidence["wide_peak_period_days"], 9.6, delta=0.25)
        self.assertGreater(evidence["power_at_twice_measured_period"], 0.9)
        self.assertLess(evidence["power_at_measured_period"], 0.15)


if __name__ == "__main__":
    unittest.main()
