"""Offline invariants for the independent 2025 MAVEN particle-event scan."""
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

try:
    import numpy as np
    import pandas as pd
    from examples.maven_sep_final_season import (
        ELECTRON, INPUT_SHA256, ION, event_orbit_counts, find_events, run, sep_metrics, summarize_hours,
    )
    HAS_ANALYSIS_EXTRA = True
except ImportError:
    HAS_ANALYSIS_EXTRA = False


@unittest.skipUnless(HAS_ANALYSIS_EXTRA, 'needs the optional numpy+pandas analysis environment')
class MavenSepFinalSeasonTests(unittest.TestCase):
    def make_data(self, count=400):
        times = pd.date_range("2025-08-25T20:00:00Z", periods=count, freq="8s")
        df = pd.DataFrame({"time": times})
        for channel in ION + ELECTRON:
            df[channel] = 2000.0 if channel in ION else 5000.0
            # This is a 1-sigma uncertainty, not a 0=good quality flag.
            df[channel + "_quality"] = 15.0
        return df

    def test_both_species_and_three_views_required(self):
        df = self.make_data(4)
        result = sep_metrics(df)
        self.assertTrue(result["joint"].all())
        self.assertTrue(np.allclose(result["ion_median"], 2000.0))
        # A single hot ion channel is never enough.
        for channel in ION[1:]:
            df.loc[0, channel] = np.nan
        result = sep_metrics(df)
        self.assertFalse(result.loc[0, "joint"])
        # An electron low-SNR result rules out a joint detection.
        df = self.make_data(4)
        for channel in ELECTRON:
            df.loc[0, channel + "_quality"] = 10000
        self.assertFalse(sep_metrics(df).loc[0, "joint"])

    def test_invalid_and_fills_dont_become_events(self):
        df = self.make_data(3)
        df.loc[0, ION[0]] = -1e31
        df.loc[1, ELECTRON[0]] = np.inf
        df.loc[2, ELECTRON[0] + "_quality"] = 0
        out = sep_metrics(df)
        self.assertTrue(out.loc[0, "joint"])  # other 3 views valid
        self.assertTrue(out.loc[1, "joint"])
        self.assertTrue(out.loc[2, "joint"])  # other 3 SNR views

    def test_sustained_two_hour_orbit_control(self):
        df = self.make_data(1350)  # 3 UTC hours with 8s samples
        sig = sep_metrics(df)
        sig["highalt"] = True
        sig["joint_highalt"] = sig["joint"]
        sig["spice_orbit_number"] = 23456
        h = summarize_hours(sig)
        self.assertEqual(int(h["event_hour"].sum()), 3)
        events = find_events(h)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["selected_hours"], 3)
        sig["highalt"] = False
        sig["joint_highalt"] = False
        self.assertEqual(int(summarize_hours(sig)["event_hour"].sum()), 0)

    def test_orbit_seen_only_in_gap_is_excluded_from_selected_hours(self):
        hours = pd.date_range("2025-08-23", periods=3, freq="h", tz="UTC")
        joint_orbits = pd.DataFrame({"hour": hours, "spice_orbit_number": [100, 101, 102]})
        selected = pd.Series([hours[0], hours[2]])
        counts = event_orbit_counts(joint_orbits, selected, hours[0], hours[2])
        self.assertEqual(counts, {
            "distinct_joint_orbits": 2,
            "full_window_distinct_joint_orbits": 3,
            "full_window_hours": 3,
        })

    def test_substituted_input_stops_replay_before_reading_or_writing_results(self):
        for instrument in tuple(INPUT_SHA256):
            with self.subTest(instrument=instrument), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                output, daily = root / "summary.json", root / "daily.csv"
                # Exercise the real downloader and digest check for each frozen instrument.
                with patch.dict(INPUT_SHA256, {instrument: INPUT_SHA256[instrument]}, clear=True), \
                     patch("urllib.request.urlopen", return_value=io.BytesIO(b"changed input")), \
                     patch("pandas.read_parquet") as read_parquet:
                    with self.assertRaisesRegex(ValueError, f"{instrument} SHA-256 mismatch"):
                        run(output, daily)
                read_parquet.assert_not_called()
                self.assertFalse(output.exists())
                self.assertFalse(daily.exists())

    def test_isolated_hour_rejected(self):
        df = self.make_data(200)
        sig = sep_metrics(df)
        sig["highalt"] = True
        sig["joint_highalt"] = True
        sig["spice_orbit_number"] = 23456
        self.assertEqual(find_events(summarize_hours(sig)), [])


if __name__ == "__main__":
    unittest.main()
