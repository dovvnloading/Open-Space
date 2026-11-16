"""Reproduce independent Mars DONKI forecast arithmetic without external services."""
import unittest
from datetime import timedelta

from examples.maven_cme_forecast_alignment import CME_FORECASTS, check, utc


class PreEventMavenForecastChecks(unittest.TestCase):
    def sample(self):
        return {
            "top_observed_hourly_fluxes": [
                {"hour_utc": "2025-08-24T05:00:00+00:00", "ion_median": 157000, "electron_median": 616500},
                {"hour_utc": "2025-08-24T04:00:00+00:00", "ion_median": 112000, "electron_median": 244500},
            ],
            "ranked_events": [
                {"start_utc": "2025-08-22T21:00:00+00:00"}
            ],
        }

    def test_original_mars_forecasts_precede_observed_peak(self):
        forecast = check(self.sample())
        self.assertTrue(forecast["august_peak_within_prior_cme21_forecast_window"])
        self.assertEqual(len(forecast["all_original_pre_event_mars_forecasts"]), 5)
        self.assertEqual(
            forecast["august_21_cme_forecast_window_utc"],
            ["2025-08-24T00:21:00+00:00", "2025-08-24T09:40:00+00:00"],
        )
        self.assertEqual(
            [x["cme_id"] for x in CME_FORECASTS].count("2025-08-20T05:12:00-CME-001"), 1
        )

    def test_05_00_vs_04_20_is_exactly_plus_40_minutes(self):
        f = check(self.sample())["all_original_pre_event_mars_forecasts"]
        self.assertAlmostEqual(
            next(x["hours_from_predicted_mars_arrival_to_peak"] for x in f
                 if x["forecast_mars_utc"] == "2025-08-24T04:20:00Z"),
            2 / 3,
            places=3,
        )
        self.assertAlmostEqual(
            next(x["hours_from_predicted_mars_arrival_to_peak"] for x in f
                 if x["forecast_mars_utc"] == "2025-08-23T22:00:00Z"),
            7.0,
        )

    def test_later_observed_peak_does_not_match_prior_window(self):
        s = self.sample()
        s["top_observed_hourly_fluxes"][0]["hour_utc"] = "2025-08-30T05:00:00Z"
        f = check(s)
        self.assertFalse(f["august_peak_within_prior_cme21_forecast_window"])

    def test_source_dates_originally_submitted_before_peak(self):
        peak = utc(self.sample()["top_observed_hourly_fluxes"][0]["hour_utc"])
        self.assertTrue(all(utc(x["analysis_submitted_utc"]) < peak for x in CME_FORECASTS))
        self.assertTrue(all(utc(x["cme_detected_utc"]) <= utc(x["analysis_submitted_utc"]) for x in CME_FORECASTS))


if __name__ == "__main__":
    unittest.main()
