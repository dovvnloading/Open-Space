"""Run the synthetic calibration and geometry regression with unittest."""
import unittest

import numpy as np

from wake_filter import compatible, detect
from wake_quality import extended_profile
from wake_search import add_wake


class WakeTests(unittest.TestCase):
    def test_blind_lines_and_negative_control(self):
        rng = np.random.default_rng(17)
        noise = rng.normal(size=(160,160)).astype("f4")
        valid = np.ones_like(noise,dtype=bool)
        for angle in (0.,33.,92.,163.):
            truth = dict(x=79.3,y=77.7,pa_image_deg=angle,length_px=60.,fwhm_px=2.)
            image = add_wake(noise,truth,2.)
            rows = detect(image,valid,[60.],[2.],.1)
            matches = [r for r in rows if compatible(truth,r)]
            self.assertTrue(matches, (angle,rows))
            self.assertTrue(any(extended_profile(image,valid,r)["profile_gate"] for r in matches))
        self.assertFalse(detect(noise,valid,[60],[2],.1))
        self.assertFalse(detect(image,~valid,[60],[2],.1))

    def test_surface_brightness_and_length_normalization(self):
        image = np.zeros((160,160),dtype="f4")
        length, width = 60., 3.
        # Integral of the Gaussian over the FWHM strip is erf(sqrt(ln(2))) of total light.
        mean = .005
        from math import erf, log, sqrt
        total = mean*length*width/erf(sqrt(log(2)))
        for angle in (0,33,145):
            truth = dict(x=79.13,y=78.33,pa_image_deg=angle,length_px=length,fwhm_px=width)
            actual = add_wake(image,truth,mean).sum()
            self.assertAlmostEqual(float(actual/total),1,places=2)

    def test_two_distinct_parallel_segments(self):
        a = dict(x=100,y=100,pa_image_deg=179,length_px=80,fwhm_px=2)
        self.assertTrue(compatible(a, dict(a,pa_image_deg=0)))
        self.assertFalse(compatible(a, dict(a,y=110)))
        self.assertFalse(compatible(a, dict(a,x=150)))
        self.assertFalse(compatible(a, dict(a,pa_image_deg=90)))


if __name__ == "__main__":
    unittest.main()
