"""Offline synthetic checks for the new GWTC-lensing search, independent of Zenodo."""
import json
import math
import pathlib
import sys
import unittest
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import run


class ScientificChecks(unittest.TestCase):
    def test_prior_finite(self):
        x=np.array([[math.log(30)], [0], [0], [0], [0], [0]], dtype=float)
        self.assertTrue(np.isfinite(run.log_prior(x)).all())

    def test_sis_delay_is_physical(self):
        t,w=run.sis_delays([.3, .8, 1.2], num=5000)
        self.assertGreater(np.mean(t > 0), .999)
        self.assertGreater(np.mean((t>.01)&(t<365)), .5)
        self.assertTrue(np.isfinite(w).all())
        self.assertTrue((w >= 0).all())

    def test_same_sky_more_consistent(self):
        try:
            import healpy
        except ImportError:
            self.skipTest("healpy available in runner, not offline")
        rng=np.random.default_rng(42)
        n=1000
        a={"ra":rng.normal(1., .03, n), "dec":rng.normal(.4,.03,n)}
        b={"ra":rng.normal(1., .03, n), "dec":rng.normal(.4,.03,n)}
        away={"ra":rng.normal(4., .03,n), "dec":rng.normal(-.7,.03,n)}
        self.assertGreater(run.sky_bayes(a,b),run.sky_bayes(a,away)+5)

    def test_detector_frame_chirp_proxy(self):
        e={"params":{"chirp_mass_source":{"best":20., "upper_error":2.0, "lower_error":-2.0},
                     "redshift":{"best":.5,"upper_error":.2,"lower_error":-.2},
                     "chi_eff":{"best":.1,"upper_error":.3,"lower_error":-.3},
                     "p_astro":{"best":.99},"far":{"best":.01},
                     "mass_2_source":{"best":15}},
           "catalog":"GWTC-4.0","gps":1380000000}
        files={"files":{"GWTC-4.0":{"GW230829_000000":{"key":"file.hdf5"}},
                        "GWTC-2.1":{},"GWTC-3":{}}}
        pop=run.population({"GW230829_000000":e},files)
        self.assertEqual(len(pop),1)
        self.assertAlmostEqual(pop[0]["mc_det"],30.)

    def test_control_pairs_are_not_claimed(self):
        self.assertIn(("GW231113_200417","GW231223_202619"),run.CONTROL_PAIRS)
        self.assertEqual(len(set(tuple(sorted(pair)) for pair in run.CONTROL_PAIRS)),len(run.CONTROL_PAIRS))


if __name__=="__main__":
    unittest.main()
