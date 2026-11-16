"""Regression checks for delay sign, crop timing, vetoes, and shifted exposure."""
import unittest

import numpy as np

import core


class SearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fs = 2048
        cls.n = 64
        rng = np.random.default_rng(1404)
        cls.h, cls.l = rng.normal(size=(2, cls.n * cls.fs)) * 1e-22
        # Multiple frequencies resolve the sign / half-cycle lag degeneracy.
        for frequency in (116, 151, 202):
            cls.h += core.sine_gaussian(cls.fs, cls.n, 28.44, frequency, 12, 1e-22)
            cls.l -= core.sine_gaussian(cls.fs, cls.n, 28.446, frequency, 12, 9e-23)
        cls.hmask = np.ones(cls.n, dtype=bool)
        cls.lmask = np.ones(cls.n, dtype=bool)
        cls.hmask[44] = False
        cls.lmask[40] = False
        cls.hp = core.prepare(cls.h, cls.fs, cls.hmask)
        cls.lp = core.prepare(cls.l, cls.fs, cls.lmask)

    def test_exposure_and_shift_sign(self):
        slides = core.scan(self.hp, self.lp, (0, 5, -5, 65))
        for lag, slide in slides.items():
            trusted = [t for t in range(self.n) if 0 <= t + lag < self.n
                       and self.hp.live[t] and self.lp.live[t + lag]]
            self.assertEqual(slide.live_seconds, len(trusted))
            self.assertEqual(slide.nominal_tiles, 88 * len(trusted))
            for tr in slide.triggers:
                self.assertIn(int(tr.time), trusted)
        self.assertTrue(np.all(~self.hp.live[43:46]))
        self.assertTrue(np.all(~self.lp.live[39:42]))

    def test_delayed_opposite_polarity_positive_control(self):
        slide = core.scan(self.hp, self.lp)[0]
        tr = max(slide.triggers, key=lambda x: x.score)
        self.assertLess(abs(tr.time - 28.44), 0.13)
        self.assertEqual(tr.sign, -1)
        self.assertLess(abs(tr.lag_seconds - 0.006), 0.001)
        self.assertEqual(tr.band, (96, 224))
        result = core.confirm(self.h, self.l, self.fs, tr.time, band=tr.band,
                              duration=tr.duration, good_h=self.hmask, good_l=self.lmask)
        self.assertGreater(result.coherence, 0.8)
        self.assertLess(abs(result.lag_seconds - 0.006), 0.001)
        self.assertEqual(result.sign, -1)

    def test_exact_recovery_on_interior_crop(self):
        start, stop = 24, 33
        crop = slice(start*self.fs, stop*self.fs)
        pp = [core.prepare(a[crop], self.fs, mask[start:stop], reference=ref, first_second=start)
              for a, mask, ref in ((self.h, self.hmask, self.hp), (self.l, self.lmask, self.lp))]
        tr = max(core.scan(*pp)[0].triggers, key=lambda x: x.score)
        full = max(core.scan(self.hp, self.lp)[0].triggers, key=lambda x: x.score)
        self.assertEqual(tr.time, full.time)
        self.assertAlmostEqual(tr.score/full.score, 1.0, places=5)
        self.assertEqual(tr.lag_seconds, full.lag_seconds)

    def test_nonfinite_default_and_degenerate_runs(self):
        x = self.h.copy()
        x[12*self.fs] = np.nan
        result = core.prepare(x, self.fs, np.ones(self.n, dtype=bool))
        self.assertTrue(np.all(~result.live[11:14]))
        dead = core.prepare(np.zeros(32*self.fs), self.fs, np.ones(32, dtype=bool))
        self.assertFalse(dead.live.any())
        slide = core.scan(dead, dead)[0]
        self.assertEqual((slide.live_seconds, slide.triggers), (0, []))

    def test_software_amplitude_and_bad_arguments(self):
        wave = core.sine_gaussian(self.fs, 2, 0.8, 154, 8.9, 1e-22)
        self.assertAlmostEqual(np.sqrt(sum(wave**2)/self.fs)/1e-22, 1., places=12)
        with self.assertRaises(ValueError):
            core.scan(self.hp, self.lp, (0, 1))
        with self.assertRaises(ValueError):
            core.prepare(self.h, self.fs, np.ones(self.n) * 2)


if __name__ == "__main__":
    unittest.main()
