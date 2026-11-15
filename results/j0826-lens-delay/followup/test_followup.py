"""Self-contained checks for new follow-up study; never edit parent PR artifacts."""
import csv
import sys
import unittest
from collections import Counter
from math import cos, hypot, radians
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from gaussian_process import observations, DelayedGP, optimize
from free_knots import knots, objective, SOURCE_RANGE, MICRO_RANGE
from sis_sanity import sis, GMAG_A, GMAG_B


class FollowupEvidenceTests(unittest.TestCase):
    def test_sources_are_fixed_upstream_and_count_matches(self):
        d = observations("ztf",width=14)
        self.assertEqual(len(d),239)
        self.assertEqual(len(observations("gaia")),76)
        self.assertEqual(len({str(x) for x in d[:,0]}), len(d))

    def test_same_time_image_model_is_weaker_for_locked_intrinsic_gp(self):
        gp=DelayedGP(observations("ztf",width=14))
        p=optimize(gp,210).x
        self.assertLess(gp.score(210,p),gp.score(0,p))
        self.assertLess(gp.score(210,p),gp.score(365,p))

    def test_season_holdout_is_not_uniform_confirmation(self):
        gp=DelayedGP(observations("ztf",width=14))
        for split in (59600,60200):
            early=DelayedGP(gp.d[gp.t<split])
            p=optimize(early,210).x
            d210=gp.heldout_conditional_nll(split,210,p)
            d0=gp.heldout_conditional_nll(split,0,p)
            if split==59600:self.assertLess(d210,d0)
            else:self.assertGreater(d210,d0)

    def test_joint_microlensing_spline_prefers_delayed_for_fixed_knots(self):
        q,m=knots(140,SOURCE_RANGE),knots(700,MICRO_RANGE)
        self.assertLess(objective(210,q,m),objective(0,q,m))

    def test_ptf_third_survey_lacks_b_and_all_uncalibrated(self):
        with (HERE/"ptf-epoch-query.csv").open(newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),27)
        counts=Counter(r["oid"] for r in rows)
        self.assertEqual(counts,{"52342070001397":17,"52342070001419":10})
        self.assertEqual({r["photcalflag"] for r in rows},{"0"})
        A=(126.53489916666669,70.04488027777778)
        B=(126.53333583333334,70.04640916666666)
        for oid in counts:
            r=next(r for r in rows if r["oid"]==oid)
            ra,dec=float(r["ra"]),float(r["dec"])
            da=3600*hypot((ra-A[0])*cos(radians(A[1])),dec-A[1])
            db=3600*hypot((ra-B[0])*cos(radians(B[1])),dec-B[1])
            self.assertGreater(db,4.0)
            self.assertGreater(da,0.8)

    def test_toy_geometry_not_a_mass_measurement(self):
        t,beta,sigma,mass,ratio=sis(.5)
        self.assertAlmostEqual(t,1530.5,delta=1.)
        self.assertAlmostEqual(beta,.137,delta=.002)
        self.assertGreater(10**(.4*(GMAG_B-GMAG_A)),ratio*3)


if __name__=="__main__":
    unittest.main()
