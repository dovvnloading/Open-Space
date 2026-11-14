"""Numerical checks using constructed velocities, independent of any future measurement."""

import unittest

import numpy as np
import pandas as pd
from astropy.time import Time

from components import score
from orbit import CLOCK_SCALE, DOSSIER, EPOCH, NOMINAL, predict_jd

import sys

sys.path.insert(0, str(DOSSIER))
from measure import gaia_delay


class ComponentTest(unittest.TestCase):
    def test_velocity_is_derivative_of_published_light_travel_time(self):
        center = EPOCH + np.linspace(200, 500, 35)
        minus = Time(center - 0.001, format="jd", scale="tcb").tdb.jd - 2457000
        plus = Time(center + 0.001, format="jd", scale="tcb").tdb.jd - 2457000
        derivative = (
            (gaia_delay(plus) - gaia_delay(minus)) * 299792.458 / (0.002 * 86400)
        )
        np.testing.assert_allclose(
            derivative, predict_jd(center) - NOMINAL[2], atol=5e-5
        )

    def test_both_extreme_but_finite_mass_ratios(self):
        jd = EPOCH + np.arange(16) * 9.3
        x = predict_jd(jd) - NOMINAL[2]
        relative_velocity = np.sin(np.arange(len(jd)) * 1.7) * 170
        for w in [0.004, 0.998]:
            pair = 12 + CLOCK_SCALE * x
            table = pd.DataFrame(
                dict(
                    exposure_id=np.arange(len(jd)),
                    jd_tcb=jd,
                    rv_a=pair + w * relative_velocity,
                    rv_b=pair - (1 - w) * relative_velocity,
                    rv_c=13 + x,
                    err_a=2.0,
                    err_b=3.0,
                    err_c=1.0,
                )
            )
            result = score(table)["inner_pair"]["predicted"]
            self.assertAlmostEqual(result["weight_b"], w, places=6)
            self.assertFalse(result["single_line_boundary"])

    def test_unequal_mass_binary_and_reversed_velocity(self):
        jd = EPOCH + np.repeat([200, 300, 380], 4) + np.tile(np.arange(4) * 0.08, 3)
        x = predict_jd(jd) - NOMINAL[2]
        drift = np.tile([98, 12, -90, -32], 3)
        q = 0.62
        for scale in (CLOCK_SCALE, -CLOCK_SCALE, 0.0):
            a = 12 + scale * x + drift
            b = 12 + scale * x - drift / q
            table = pd.DataFrame(
                dict(
                    exposure_id=np.arange(len(jd)),
                    jd_tcb=jd,
                    rv_a=a,
                    rv_b=b,
                    rv_c=14 + x,
                    err_a=2.0,
                    err_b=4.0,
                    err_c=1.0,
                )
            )
            result = score(table)
            self.assertTrue(result["free_pair_model_identifiable"])
            inner = result["inner_pair"]["free"]
            self.assertAlmostEqual(inner["mass_ratio_b_over_a"], q, places=5)
            self.assertAlmostEqual(inner["outer_scale"], scale, places=5)
            self.assertAlmostEqual(
                result["third_line"]["free"]["outer_scale"], 1.0, places=5
            )
        with self.assertRaisesRegex(ValueError, "unique"):
            score(pd.concat([table, table.iloc[:1]]))
        with self.assertRaisesRegex(ValueError, "single predicted"):
            score(table.assign(jd_tcb=jd[0]))


if __name__ == "__main__":
    unittest.main()
