"""Regression checks for the absorption coefficient and the spectral likelihood profile."""

import unittest

import numpy as np
from scipy.optimize import minimize

from fit_pha import profile_shape, profile_stat, sigma


class SpectrumProfileTests(unittest.TestCase):
    def test_morrison_mccammon_5_kev(self):
        self.assertAlmostEqual(sigma(5.0) * 1.0e24, 3.5252, places=10)

    def test_profile_matches_independent_joint_minimization(self):
        matrix = np.eye(4)
        absorption = np.array([1.5, 0.9, 0.4, 0.1])
        continuum = np.array([10.0, 7.0, 3.0, 1.0])
        counts = np.array([3.0, 12.0, 20.0, 7.0])
        background = np.array([1.0, 3.0, 2.0, 1.0])
        alpha = 0.2

        def joint(params):
            total, column = params
            unit = np.exp(-column * absorption) * continuum
            return profile_stat(counts, background, alpha, total * unit / unit.sum())

        for bounds in [(0.0, 10.0), (0.5, 0.5)]:
            reference = minimize(
                joint,
                [40.0, 0.5],
                bounds=[(0.0, 100.0), bounds],
                method="Powell",
                options={"xtol": 1.0e-10, "ftol": 1.0e-12},
            )
            self.assertTrue(reference.success)
            for grid in (9, 18):
                with self.subTest(bounds=bounds, grid=grid):
                    value = profile_shape(
                        continuum,
                        matrix,
                        absorption,
                        counts,
                        background,
                        alpha,
                        [(-8, 16), bounds],
                        grid,
                    )
                    self.assertAlmostEqual(value, reference.fun, places=7)


if __name__ == "__main__":
    unittest.main()
