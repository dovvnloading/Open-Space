"""Deterministic engineering regressions; the astronomical controls are measured separately."""
import unittest

import numpy as np

from recovery import geometry, inject, matches, occurrence, overlap
from wd_search import search


class SearchChecks(unittest.TestCase):
    def test_overlap_and_dilution(self):
        np.testing.assert_allclose(overlap(np.array([0., 3.]), .6), [.36, 0.])
        np.testing.assert_allclose(overlap(np.array([0., 3.]), 1.5), [1., 0.])
        t = np.arange(3000)*20/86400
        visit = (t, np.ones(len(t)), np.ones(len(t)), 20., .4)
        a,k,b,w = geometry(.25, 11.2, .60, .012, 0.)
        f = inject(visit, .25, t[100], a,k,b,w)
        self.assertAlmostEqual(f[100], .60)
        self.assertAlmostEqual(max(f), 1.)

    def test_phase_wrap_in_all_cadences(self):
        period = 1.40794
        for dt in (20., 120., 200.):
            t = np.arange(0.,12.,dt/86400.) + 4000.
            t = t[(t < 4004.1) | (t > 4004.6)]
            flux = np.random.default_rng(7).normal(1., .015, len(t))
            epoch = t[0]-.001
            f = flux - .15*(abs((t-epoch+.5*period)%period-.5*period) < .003)
            true = (t, f, np.ones_like(f), dt, 1.)
            hits = search(t, f, dt)['periodic']
            self.assertTrue(any(matches(row,true,period,epoch,.006) for row in hits))
            self.assertFalse(matches(dict(period_days=.5*period,epoch_btjd=epoch),
                                     true,period,epoch,.006))

    def test_monotone_confidence_bound(self):
        s = [dict(period_days=1.,radius_earth=2.,target=str(j),trials=100,
                  recovered=95,possible_host=0,geometry_probability=.01) for j in range(20)]
        initial = occurrence(s)[0]
        self.assertLess(initial['transiting_host_upper95'], .25)
        self.assertIsNone(initial['intrinsic_host_upper95'])
        s[0]['possible_host'] = 1
        self.assertGreater(occurrence(s)[0]['transiting_host_upper95'], initial['transiting_host_upper95'])


if __name__ == '__main__':
    unittest.main()
