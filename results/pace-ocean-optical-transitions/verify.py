"""Offline replay of mission headline measurements from the committed tables."""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def table(path):
    return pd.read_csv(ROOT / path)


def check_close(actual, expected, tolerance):
    np.testing.assert_allclose(actual, expected, atol=tolerance, rtol=0)


def check_search():
    coverage = table("platforms/hyperspectral_coverage.csv")
    assert coverage.matchups.sum() == coverage.complete_oci_61.sum() == 932
    assert coverage.complete_field_61.sum() == 118
    intervals = table("platforms/fixed_platform_intervals.csv")
    assert len(intervals) == 644
    assert intervals[["SAT_blue_delta", "FIELD_blue_delta"]].dropna().shape[0] == 599
    comparisons = table("florida/regional_changes.csv")
    passing = comparisons[comparisons.coverage_pass]
    assert len(comparisons) == 684 and len(passing) == 9
    assert passing.paired_points.sum() == 256
    assert table("florida/event_spectra.csv").wavelength_nm.nunique() == 149
    matched = table("sentinel3/olci_footprints_summary.csv").query(
        "selection == 'all_positive_bands'"
    )
    assert matched.collocations.sum() == 236
    event = matched[
        (matched.region == "florida_atlantic") & (matched.start == "2025-05-13")
    ].iloc[0]
    assert event.collocations == 14
    check_close(
        [event.oci_median_delta_C, event.olci_median_delta_C],
        [-0.23659, -0.00157],
        0.00002,
    )


def check_adriatic():
    spectra = table("adriatic/spectra.csv")
    for instrument, expected in [("oci_rrs", -0.1736), ("field_rrs", -0.3367)]:
        w = spectra.pivot(index="date", columns="wavelength_nm", values=instrument)
        c = np.log(w[620]) - (45 * np.log(w[560]) + 60 * np.log(w[665])) / 105
        check_close(c.iloc[1] - c.iloc[0], expected, 0.0001)
    cells = table("sentinel3/aaot_paired_spectra.csv").query(
        "distance_km <= 50 and depth_m > 0"
    )

    def curvature(day):
        blue, center, red = [cells[f"Rrs{b}_{day}"] for b in (560, 620, 665)]
        return np.log(center) - (45 * np.log(blue) + 60 * np.log(red)) / 105

    delta = curvature("2025-11-07") - curvature("2025-11-06")
    check_close(delta, cells.delta_C, 0.000001)
    assert len(cells) == 268 and (delta < 0).sum() == 225
    check_close(delta.median(), -0.08260, 0.00001)
    step = 1 / 24
    areas = (
        6371**2
        * np.deg2rad(step)
        * (
            np.sin(np.deg2rad(cells.lat + step / 2))
            - np.sin(np.deg2rad(cells.lat - step / 2))
        )
    )
    check_close([areas.sum(), areas[delta < 0].sum()], [4050, 3400], 1)
    control = table("sentinel3/aaot_regional_seasonal.csv")
    assert len(control) == 58
    extreme = control[control.historical_median <= control.matched_event_median]
    assert set(extreme.reference_start) == {"2019-10-25", "2026-10-30"}
    context = table("sentinel3/aaot_spectral_context.csv")
    rebound = context[
        (context.start == "2025-11-05") & (context.end == "2025-11-07")
    ].iloc[0]
    assert rebound.paired_all_bands == 262
    check_close(np.expm1(rebound.median_log_change_620), 0.0059, 0.0001)
    casts = table("veneto/veneto_casts.csv")
    assert len(casts) == 29
    check_close(casts[["turbidity_FTU", "chl_a_ug_L"]].median(), [2.67, 0.56], 0.001)
    daily = table("adriatic/field_daily_event.csv").set_index("date").iloc[:, 0]
    fall = table("adriatic/field_daily_increments.csv")
    fall = fall[fall.date < "2025-01-01"]
    assert len(fall) == 151
    assert (
        fall.delta_C560_620_667 <= daily["2025-11-07"] - daily["2025-11-06"]
    ).sum() == 5


if __name__ == "__main__":
    check_search()
    check_adriatic()
    print(
        "Verified search denominators, Florida collocations, Adriatic spectra, sampled area, historical comparisons, field casts and two-day recovery."
    )
