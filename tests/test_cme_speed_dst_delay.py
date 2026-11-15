"""Synthetic data-replay controls for CME -> Dst latency inference."""
import importlib.util
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[1]/"examples"/"cme_speed_dst_delay.py"
s=importlib.util.spec_from_file_location("cme_speed_dst_delay",path)
mod=importlib.util.module_from_spec(s)
s.loader.exec_module(mod)

HEADER=("No,Peak Date/Time,DST Peak (nT),Bz Peak(nT),"
        "icme onsett (yyyy-mm-dd hh:mm),Cme speed (km/s),ACTIVE REGION")
ROWS=("1,2001-2-3 7:00,-120,-11,2001-2-3 1:30,1500,AR 999\n"
      "2,2001-2-10 12:00,-123,-11,CIR,NA,NA\n"
      "3,2001-2-12 9:00,-140,-20,2001-2-11 6:00,NA,Quiet Region\n"
      "4,2001-2-18 17:00,-110,-8,2001-2-16 15:00,400,Quiet Region\n")

class CMELatencyTests(unittest.TestCase):
    def test_variable_width_utc_dates(self):
        self.assertEqual(mod.parse_utc("2001-2-3 7:00").hour,7)
        self.assertIsNone(mod.parse_utc("CIR"))
        self.assertIsNone(mod.parse_utc("NA"))

    def test_excludes_cir_and_missing_cme_speed(self):
        events,counts=mod.parse_catalog(HEADER+"\n"+ROWS)
        self.assertEqual(counts["total_storm_rows"],4)
        self.assertEqual(counts["no_shock_timestamp"],1)
        self.assertEqual(counts["missing_cme_speed"],1)
        self.assertEqual(len(events),2)
        self.assertAlmostEqual(events[0]["delay_hours"],5.5)

    def test_strong_anticorrelation_is_not_assumed(self):
        self.assertAlmostEqual(mod.pearson([(10,40),(20,30),(30,20)]),-1)
        self.assertAlmostEqual(mod.pearson([(10,20),(20,30),(30,40)]),1)
        with self.assertRaises(ValueError):
            mod.pearson([(1,2),(1,3),(1,4)])

    def test_median_gap_direction_and_no_threshold_confusion(self):
        f=lambda speed,delay:{"cme_speed_km_s":speed,"delay_hours":delay}
        events=[f(600,30),f(400,20),f(1100,7),f(1200,9),f(800,17)]
        self.assertEqual(mod.median_gap(events),17)
        self.assertEqual(mod.cohort(events)["slow_cme_lt700_n"],2)
        self.assertEqual(mod.cohort(events)["fast_cme_ge1000_n"],2)

    def test_separate_shock_episodes(self):
        d=[{"shock_epoch":0},{"shock_epoch":10*3600},
           {"shock_epoch":70*3600},{"shock_epoch":73*3600}]
        e=mod.independent_shock_episodes(d)
        self.assertEqual([x["shock_epoch"] for x in e],[0,73*3600])

    def test_fixed_permutation_seed(self):
        d=[{"year":1997+(i//4),"cme_speed_km_s":300+150*i,
            "delay_hours":40-i*i/10} for i in range(8)]
        a=mod.era_permutation(d,60,32)
        b=mod.era_permutation(d,60,32)
        self.assertEqual(a,b)
        self.assertGreaterEqual(a["two_sided_empirical_p"],round(1/61,7))

if __name__=="__main__":
    unittest.main()
