"""Known-answer scientific controls for local fixed-ephemeris measurements.

These fixtures are explicitly synthetic and cannot enter the observed register.
"""
import importlib.util
from pathlib import Path
import unittest

try:
    import numpy as np
except ImportError:
    np = None

PATH = Path(__file__).resolve().parents[1] / 'results/exoplanet-search-20261125/vet.py'
if np is not None:
    SPEC = importlib.util.spec_from_file_location('exoplanet_vet', PATH)
    vet = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(vet)


@unittest.skipIf(np is None, 'NumPy is an optional research dependency')
class LocalTransitControls(unittest.TestCase):
    @staticmethod
    def curve(depth=.002, alternating=False):
        rng = np.random.default_rng(125269)
        time = np.arange(0, 24, .002)
        period, epoch, duration = 2.5, 1., .10
        cycle = np.rint((time - epoch) / period).astype(int)
        in_box = np.abs(time - epoch - cycle * period) < duration / 2
        baseline = 1 + .001 * time
        depth_array = depth * (1 + (cycle % 2) if alternating else np.ones(len(time)))
        flux = baseline * (1 - depth_array * in_box) + rng.normal(0, .00015, len(time))
        data = dict(time=time, pdcsap=flux, sap=flux.copy(), quality=np.zeros(len(time), int),
                    pdcsap_err=np.full(len(time), .00015), sap_err=np.full(len(time), .00015),
                    bkg=np.ones(len(time)), cx=np.ones(len(time)), cy=np.ones(len(time)))
        return data, period, epoch, duration

    def test_recovers_known_depth_on_sloping_unfiltered_flux(self):
        data, p, t0, d = self.curve()
        events = vet.all_events(data, p, t0, d)
        result = vet.aggregate(events)
        self.assertGreaterEqual(result['covered_events'], 8)
        self.assertAlmostEqual(result['equal_event_depth'], 2000, delta=70)
        self.assertGreater(result['drop_one_min_equal_event_snr'], 20)
        self.assertEqual(result['units'], 'ppm')
        # A wrong phase is independently specified and contains no injected dip.
        wrong = vet.aggregate(vet.all_events(data, p, t0 + .37 * p, d))
        self.assertLess(abs(wrong['equal_event_depth']), 80)

    def test_alternating_eclipses_are_exposed(self):
        data, p, t0, d = self.curve(alternating=True)
        result = vet.odd_even(vet.all_events(data, p, t0, d))
        self.assertAlmostEqual(result['even_minus_odd_depth'], -2000, delta=90)
        self.assertLess(result['difference_snr'], -15)

    def test_missing_flank_is_not_accepted_as_a_dip(self):
        data, p, t0, d = self.curve()
        data['quality'][(data['time'] > t0 - .4) & (data['time'] < t0 - .07)] = 1
        result = vet.measure_event(data, t0, d, p)
        self.assertFalse(result['eligible'])
        self.assertEqual(result['counts']['left'], 0)
        self.assertNotIn('depth', result)

    def test_brightenings_keep_their_negative_sign(self):
        data, p, t0, d = self.curve(depth=-.002)
        result = vet.aggregate(vet.all_events(data, p, t0, d))
        self.assertAlmostEqual(result['equal_event_depth'], -2000, delta=70)
        self.assertEqual(result['positive_events'], 0)


if __name__ == '__main__':
    unittest.main()
