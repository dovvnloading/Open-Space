import csv
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from examples import template_residuals
from open_space.dataio import digest
from open_space.holdout import screen, window

try:
    from astropy.io import fits
except ImportError:
    fits = None


@unittest.skipUnless(fits, "Astropy is required for the model-COADD integration test")
class TemplateResidualTests(unittest.TestCase):
    def test_template_polarity_threshold_coverage_and_hash_check(self):
        # An absorption trough has a strong, five-pixel Gaussian core. The fit
        # can explain all, none or a fraction. One missing model core removes
        # coverage without turning the affected observation into a residual miss.
        protocol = {"predictions": [{"id": "X", "direction": "absorption"}],
                    "velocity_half_window_kms": 150}
        evidence = {"observations": [], "protocol_sha256": "test"}
        profiles = [1.0, 0.0, 0.9, -1.0, None]
        with tempfile.TemporaryDirectory() as temp:
            inputs = Path(temp)
            (inputs / "raw").mkdir()
            for i, fraction in enumerate(profiles):
                obs_id = f"fixture-{i}"
                wave = [5000.0 + k for k in range(-40, 41)]
                line = [-12 * math.exp(-0.5 * ((w - 5000) / 1.3)**2) for w in wave]
                flux = [30.0 + v for v in line]
                model = [30.0 + (fraction or 0) * v for v in line]
                if fraction is None:
                    model[40] = float("nan")
                table = fits.BinTableHDU.from_columns(
                    [fits.Column(name="model", format="D", array=model)], name="COADD")
                raw = inputs / "raw" / (obs_id + ".fits")
                fits.HDUList([fits.PrimaryHDU(), table]).writeto(raw)
                rows = [{"wavelength": w, "flux": f, "ivar": 1.0, "mask": 0}
                        for w, f in zip(wave, flux)]
                encoded = io.StringIO()
                writer = csv.DictWriter(encoded, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
                csv_bytes = encoded.getvalue().encode()
                csv_file = obs_id + ".csv"
                (inputs / csv_file).write_bytes(csv_bytes)
                derivatives = {"source_sha256": digest(raw.read_bytes()), "csv_sha256": digest(csv_bytes),
                               "redshift": 0, "observation_id": obs_id}
                metadata = {**derivatives, "schema_version": 1, "pixels": len(rows),
                            "wavelength_frame": "observed vacuum", "wavelength_unit": "Angstrom",
                            "csv_file": csv_file, "source_url": "synthetic fixture"}
                (inputs / (obs_id + ".json")).write_text(json.dumps(metadata))
                hit = window(*screen(rows, 0), 5000., "absorption", 150)
                self.assertEqual(hit["outcome"], "hit")
                evidence["observations"].append({**derivatives, "measurements": [{"id":"X", "primary":hit}]})

            comparison = template_residuals.compare(evidence, protocol, inputs)
            records = comparison["measurements"]["X"]
            self.assertEqual([r["residual_search"] for r in records],
                             ["miss", "hit", "miss", "hit", "unsearchable"])
            for index, fraction in enumerate(profiles[:-1]):
                self.assertAlmostEqual(records[index]["model_to_observed_ratio"], fraction)
            summary = comparison["summary"][0]
            self.assertEqual(summary["observed_hits"], 5)
            self.assertEqual(summary["usable_models"], 4)
            self.assertEqual(summary["same_sign_in_model"], 2)
            self.assertEqual(summary["residual_window_covered"], 4)
            self.assertEqual(summary["residual_detections"], 2)
            with raw.open("ab") as stream:
                stream.write(b"mismatched raw data")
            with self.assertRaisesRegex(ValueError, "Input differs"):
                template_residuals.compare(evidence, protocol, inputs)
