#!/usr/bin/env python3
"""Reproduce a deliberately conservative, frozen Mars-specific forecast check.

The NASA M2M DONKI predictions below were archived in 2025, before the
measured 2025 Aug 24 SEP peak. This script measures *temporal concordance*,
not the causation, identity or significance of any heliospheric event.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

# NASA CCMC M2M DONKI CME analysis, version 1, independent WSA-ENLIL runs.
# Times/IDs can be checked manually against the cited original NASA pages.
CME_FORECASTS = (
    {
        "cme_id": "2025-08-20T05:12:00-CME-001",
        "source_url": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40665/1",
        "cme_detected_utc": "2025-08-20T05:12:00Z",
        "analysis_submitted_utc": "2025-08-20T11:51:00Z",
        "forecast_mars_utc": "2025-08-23T22:00:00Z",
        "model": "NASA WSA-ENLIL+Cone, 2.0 AU, LE",
        "type": "leading_edge",
    },
    {
        "cme_id": "2025-08-21T08:23:00-CME-001",
        "source_url": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40685/1",
        "cme_detected_utc": "2025-08-21T08:23:00Z",
        "analysis_submitted_utc": "2025-08-21T13:15:00Z",
        "forecast_mars_utc": "2025-08-24T00:21:00Z",
        "model": "NASA WSA-ENLIL+Cone, 2.0 AU, LE",
        "type": "leading_edge",
    },
    {
        "cme_id": "2025-08-21T08:23:00-CME-001",
        "source_url": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40685/1",
        "cme_detected_utc": "2025-08-21T08:23:00Z",
        "analysis_submitted_utc": "2025-08-21T13:15:00Z",
        "forecast_mars_utc": "2025-08-24T04:20:00Z",
        "model": "NASA WSA-ENLIL+Cone, 5.5 AU, LE",
        "type": "leading_edge",
    },
    {
        "cme_id": "2025-08-21T08:23:00-CME-001",
        "source_url": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40685/1",
        "cme_detected_utc": "2025-08-21T08:23:00Z",
        "analysis_submitted_utc": "2025-08-21T13:15:00Z",
        "forecast_mars_utc": "2025-08-24T09:38:00Z",
        "model": "NASA WSA-ENLIL+Cone, 2.0 AU, SH",
        "type": "shock_front",
    },
    {
        "cme_id": "2025-08-21T08:23:00-CME-001",
        "source_url": "https://kauai.ccmc.gsfc.nasa.gov/DONKI/view/CME/40685/1",
        "cme_detected_utc": "2025-08-21T08:23:00Z",
        "analysis_submitted_utc": "2025-08-21T13:15:00Z",
        "forecast_mars_utc": "2025-08-24T09:40:00Z",
        "model": "NASA WSA-ENLIL+Cone, 5.5 AU, SH",
        "type": "shock_front",
    },
)


def utc(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(timezone.utc)


def check(data: dict) -> dict:
    hours = data["top_observed_hourly_fluxes"]
    peak = max(hours, key=lambda x: x["ion_median"])
    dt_peak = utc(peak["hour_utc"])
    earliest_event = utc(data["ranked_events"][0]["start_utc"])
    cme21 = [
        x for x in CME_FORECASTS
        if x["cme_id"] == "2025-08-21T08:23:00-CME-001"
    ]
    start = min(utc(x["forecast_mars_utc"]) for x in cme21)
    end = max(utc(x["forecast_mars_utc"]) for x in cme21)
    scored = []
    for x in CME_FORECASTS:
        submitted = utc(x["analysis_submitted_utc"])
        forecast = utc(x["forecast_mars_utc"])
        if submitted >= forecast or submitted >= dt_peak:
            raise ValueError("Not a pre-event forecast: " + str(x))
        scored.append({
            **x,
            "hours_from_predicted_mars_arrival_to_peak": round(
                (dt_peak - forecast).total_seconds() / 3600.0, 4
            ),
            "hours_from_forecast_submission_to_peak": round(
                (dt_peak - submitted).total_seconds() / 3600.0, 4
            ),
        })
    return {
        "classification": "independent pre-event Mars-focused forecast temporal concordance; not shock detection or unique CME association",
        "measured_peak_hour_utc": dt_peak.isoformat(),
        "measured_peak_ion_flux": peak["ion_median"],
        "measured_peak_electron_flux": peak["electron_median"],
        "august_event_earliest_selected_hour_utc": earliest_event.isoformat(),
        "august_21_cme_forecast_window_utc": [start.isoformat(), end.isoformat()],
        "august_peak_within_prior_cme21_forecast_window": start <= dt_peak <= end,
        "competing_august_20_cme_pre_event_forecast": True,
        "all_original_pre_event_mars_forecasts": scored,
        "inference_limits": [
            "All predictions refer to arrival of modeled CME/shock structures, not necessarily onset or maximum of solar energetic particles.",
            "Both the 2025-08-20 and 2025-08-21 CMEs were separately forecast to reach Mars; peak concordance cannot identify either uniquely.",
            "These predictions were individually pre-event, but the choice to test their link to the strongest observed event was exploratory; no calibrated global coincidence p-value.",
            "Attribution requires an independent calibrated MAVEN SEP Level-2 replay, upstream MAG/SWIA shock evidence and solar particle velocity/trajectory modeling.",
            "Other Mars-directed forecasts and false-positive days were not exhaustively enumerated.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path,
        default=Path("results/maven-2025-sep-events/summary.json"),
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/maven-2025-sep-events/forecast_alignment.json"),
    )
    args = parser.parse_args()
    evidence = check(json.loads(args.input.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "peak": evidence["measured_peak_hour_utc"],
        "in_forecast_window": evidence["august_peak_within_prior_cme21_forecast_window"],
        "prediction_deltas": [(x["model"], x["hours_from_predicted_mars_arrival_to_peak"]) for x in evidence["all_original_pre_event_mars_forecasts"]]
    }, indent=2))


if __name__ == "__main__":
    main()
