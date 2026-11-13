import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

try:
    import numpy as np
except ImportError:
    np = None


@unittest.skipIf(np is None, "TESS imaging needs NumPy from the archive extra")
class ImageFollowupTests(unittest.TestCase):
    def test_known_displaced_source_on_drifting_background(self):
        from open_space import pixels
        rng = np.random.default_rng(7)
        time = np.arange(0.0, 12., .008)
        y, x = np.indices((7, 8))
        profile = np.exp(-((x-5.1)**2 + (y-3.2)**2)/(2*.62**2))
        phase = abs((time - .48 + .5) % 1. - .5)
        flux = 100 + .3*time[:, None, None] + rng.normal(0, 1, (len(time), 7, 8))
        flux -= (phase[:, None, None] < .04) * 15 * profile
        args = pixels.event_images(time, flux, np.ones_like(flux), .48, 1., .1, .008)
        report, images = pixels.summarize(*args, (2., 3.), np.ones((7, 8), bool), True)
        self.assertEqual(report["disposition"], "offset_source_followup")
        self.assertLess(np.linalg.norm(np.subtract(report["centroid_xy"], (5.1, 3.2))), .2)
        self.assertGreaterEqual(report["accepted_events"], 10)
        self.assertAlmostEqual(images["difference"][3,5], 15 * profile[3,5], delta=.5)
        self.assertIsNotNone(report["event_bootstrap_radius95_pixels"])

    def test_missing_bracket_and_zero_events(self):
        from open_space import pixels
        times = np.arange(0, 3, .01)
        times = times[(times % 1.0) < .52]
        self.assertEqual(pixels._windows(times, .48, 1., .1, .01), [])
        result, images = pixels.summarize([], [], [], [], (2, 3), np.ones((2, 3), bool))
        self.assertEqual(result["disposition"], "insufficient_events")
        self.assertEqual(images, {})

    def test_stellar_extrema_in_fixed_disjoint_cycles(self):
        from open_space import pixels
        t = np.arange(0, 13., .01)
        epoch = .5
        y, x = np.indices((7, 7))
        profile = np.exp(-((x-4.9)**2 + (y-3.2)**2)/(2*.65**2))
        flux = 100 + .1*t[:,None,None] - 8*np.cos(2*np.pi*(t-epoch))[:,None,None]*profile
        args = pixels.event_images(t, flux, np.ones_like(flux), epoch, 1., 1., .01,
                                   pixels.VARIABLE_WINDOWS)
        out, _ = pixels.summarize(*args, (2., 3.), np.ones_like(profile, dtype=bool), True)
        self.assertEqual(out["disposition"], "offset_source_followup")
        self.assertLess(np.linalg.norm(np.subtract(out["centroid_xy"], (4.9, 3.2))), .2)
        self.assertGreater(out["accepted_events"], 10)
        windows = pixels._windows(t,epoch,1.,1.,.01,pixels.VARIABLE_WINDOWS)
        indices = [i for _,_,slices in windows for sl in slices for i in range(sl.start,sl.stop)]
        self.assertEqual(len(indices), len(set(indices)))

    def test_weighted_images_require_real_pixels(self):
        from open_space import pixels
        stack = np.ones((5, 3, 3))
        stack[:4,1,1] = np.nan
        mean, error = pixels._combine(stack, np.ones_like(stack))
        self.assertTrue(np.isnan(mean[1,1]))
        self.assertTrue(np.isnan(error[1,1]))
        self.assertAlmostEqual(error[0,0], 1 / np.sqrt(5))

    def test_partner_identity_and_pinned_cache(self):
        from open_space import pixels
        root = Path(__file__).resolve().parents[1]
        task = json.loads((root/"data/pixels-wasp18.json").read_text())[0]
        name = pixels.validate({**task, "product_name": task["product_name"].replace("_tp", "_lc")})[0]
        self.assertEqual(task["product_name"], name)
        for change in [{"sector": 99}, {"tic_id": 8}, {"duration_days": .4},
                       {"product_name":"../../x.fits"}, {"period_days":float("nan")}]:
            with self.assertRaises(ValueError):
                pixels.validate({**task, **change})
        with TemporaryDirectory() as td:
            (Path(td)/name).write_bytes(b"truncated")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                pixels._acquire(Path(td), name, task["source_sha256"])
