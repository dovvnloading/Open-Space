"""Synthetic timing, coverage and target-label leakage controls."""
import importlib.util
from datetime import datetime,timedelta,timezone
from pathlib import Path
import unittest

MODULE=Path(__file__).resolve().parents[1]/"examples"/"forbush_dst_independent.py"
spec=importlib.util.spec_from_file_location("forbush_dst_independent",MODULE)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
START=datetime(2007,1,1,tzinfo=timezone.utc)

class ForbushIndependentTests(unittest.TestCase):
    def test_future_only(self):
        self.assertIsNone(mod.lag_to_next(100.01,[100,200]))
        self.assertAlmostEqual(mod.lag_to_next(100.01,[100,120]),19.99)
        self.assertIsNone(mod.lag_to_next(100,[148]))
        self.assertEqual(mod.lag_to_next(100,[147]),47)

    def test_missing_history_and_storm_veto(self):
        dst={i:0.0 for i in range(200)}
        dst[75]=None
        dst[120]=-120.0
        dst[125]=-180.0
        onset,eligible=mod.severe_onsets_and_eligibility(dst)
        self.assertEqual(onset,[]) # incomplete quiet-history coverage at 120
        dst[75]=0.0
        onset,eligible=mod.severe_onsets_and_eligibility(dst)
        self.assertEqual(onset,[120])
        self.assertNotIn(120,eligible)
        self.assertNotIn(125,eligible)
        self.assertNotIn(193,eligible)
        self.assertIn(198,eligible)

    def test_cluster_uses_first_event_before_labels(self):
        a=[{"datetime_utc":START,"onset_type":9},
           {"datetime_utc":START+timedelta(hours=71),"onset_type":1},
           {"datetime_utc":START+timedelta(hours=142),"onset_type":1},
           {"datetime_utc":START+timedelta(hours=215),"onset_type":2}]
        c=mod.cluster_events(a)
        self.assertEqual(len(c),2)
        self.assertEqual(c[0]["size"],3)
        self.assertEqual(c[0]["onset_type"],9)
        self.assertEqual(c[1]["onset_type"],2)

    def test_no_leak_from_feid_embedded_dst_or_postevent_amplitude(self):
        dst=[{"datetime":START+timedelta(hours=i),
              "dst_index_nt":(-120 if i==90 else 0)} for i in range(180)]
        event=[{"datetime_utc":START+timedelta(hours=81),
                "onset_type":9,"dst_min_nt":-5000,"fd_magnitude_pct":99}]
        a=mod.analyze(event,dst,draws=20)
        self.assertEqual(a["cohorts"]["holdout_2001_2016"]["hits_48h"],1)
        event[0]["dst_min_nt"]=5000
        event[0]["fd_magnitude_pct"]=0
        self.assertEqual(a,mod.analyze(event,dst,draws=20))

    def test_current_dst_hour_is_not_lookahead_for_eligibility(self):
        # Dst hour 90 contains a storm, but an event at 90:15 may not
        # use that hour-average to decide whether it was at-risk.
        dst=[{"datetime":START+timedelta(hours=i),
              "dst_index_nt":(-120 if i==90 else 0)} for i in range(180)]
        e={"datetime_utc":START+timedelta(hours=90,minutes=15),
           "onset_type":9}
        result=mod.analyze([e],dst,draws=20)
        self.assertEqual(result["cohorts"]["holdout_2001_2016"]["eligible_episodes"],1)
        self.assertEqual(result["cohorts"]["holdout_2001_2016"]["hits_48h"],0)

    def test_shifted_controls_preserve_month_and_eligibility(self):
        h=mod.utc_hour(START)
        e=[{"time":START+timedelta(hours=25,minutes=30)}]
        slots={h+i for i in range(744)}
        r=mod.time_shift_controls(e,slots,[h+25+7*24+12])
        self.assertEqual(r["7"]["eligible_shifted_episodes"],1)
        self.assertEqual(r["7"]["hits_48h"],1)
        self.assertEqual(r["-14"]["eligible_shifted_episodes"],0)

    def test_stratified_null_preserves_event_minute(self):
        # A hazard at hour 12 must NOT be counted for a sampled timestamp
        # at hour 12:30 (past the nominal onset).
        h=mod.utc_hour(START)
        events=[{"time":START+timedelta(hours=12,minutes=30)}]
        slots={mod.hour_to_month(h+12):[h+12]}
        r=mod.stratified_null(events,[h+12],slots,repetitions=30)
        self.assertEqual(r["mean_hits"],0)

    def test_month_matched_control_deterministic(self):
        event=[{"time":START+timedelta(hours=1)}]
        anchor=mod.utc_hour(START)
        panel={mod.hour_to_month(anchor):[anchor+x for x in (1,2,200)]}
        a=mod.stratified_null(event,[anchor+203],panel,repetitions=40,seed=4)
        self.assertEqual(a,mod.stratified_null(event,[anchor+203],panel,repetitions=40,seed=4))
        self.assertEqual(a["draws"],40)

    def test_sentinel_and_duplicate_hour(self):
        rows=[{"datetime":START,"dst_index_nt":99999},
              {"datetime":START+timedelta(hours=1),"dst_index_nt":-100}]
        self.assertEqual(mod.make_hourly_dst(rows)[mod.utc_hour(START)],None)
        self.assertEqual(mod.make_hourly_dst(rows)[mod.utc_hour(START)+1],-100)
        with self.assertRaises(ValueError):
            mod.make_hourly_dst(rows+[rows[0]])

if __name__=="__main__":
    unittest.main()
