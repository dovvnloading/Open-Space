import math
import unittest

from publication import DURATIONS, instrument_rate_bounds


class ConservativeDetectedRate(unittest.TestCase):
    def test_zero_and_positive_counts_share_exact_poisson_construction(self):
        exposure = {str(d): 2e6 for d in DURATIONS}
        exposure[str(DURATIONS[-1])] = 0.
        empty = instrument_rate_bounds(0, exposure)
        self.assertAlmostEqual(empty[str(DURATIONS[0])], -math.log(.05)/2)
        self.assertIsNone(empty[str(DURATIONS[-1])])
        one = instrument_rate_bounds(1, exposure)[str(DURATIONS[0])] * 2
        self.assertAlmostEqual(math.exp(-one) * (1+one), .05)
        with self.assertRaises(ValueError):
            instrument_rate_bounds(-1, exposure)
