"""Additional clock tests using frozen PR #23 inputs and a distinct 2024 burst table.

Run: python results/nicer-magnetar-microbursts/timescale-forensics/verify.py
No third-party packages or downloads are needed. The clock conversion for the
local predicted times follows NICER's documented post-2017 calibrated TT epoch.
It is not a remeasurement of source photons or a new Fermi test.
"""
import csv
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
WINDOWS = HERE.parent / "catalog-clock" / "windows.csv"
PAIRS = HERE / "shao_2024_table4_pairs.csv"
UTC_EPOCH = datetime.fromisoformat("2014-01-01T00:00:00")
TT_MINUS_UTC_2017_PLUS = 69.184
EPOCH_TT_MINUS_UTC_2014 = 67.184
LEAPINIT = 2.0


def rows(path):
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def test_time_references():
    # NASA NICER timing documentation; 2014 UTC epoch expressed as TT.
    mjdreff = 0.000777592592592593
    assert abs(mjdreff * 86400.0 - EPOCH_TT_MINUS_UTC_2014) < 1e-9
    assert abs(EPOCH_TT_MINUS_UTC_2014 + LEAPINIT
               - TT_MINUS_UTC_2017_PLUS) < 1e-12
    assert abs(32.184 + 37.0 - TT_MINUS_UTC_2017_PLUS) < 1e-12
    assert abs(32.184 + 35.0 - EPOCH_TT_MINUS_UTC_2014) < 1e-12
    print("NASA header-model checks: 67.184-s reference + 2 leap seconds = 69.184 s")


def year_strata(entries):
    groups = defaultdict(lambda: dict(n=0, wins=0, losses=0, ties=0,
                                      count_at_label=0, count_at_plus69=0,
                                      obs=set(), sources=set()))
    for row in entries:
        if not row["narrow_literal"] or not row["narrow_plus69"]:
            continue  # Require complete paired exposure, as in the source file.
        year = row["start"][:4]
        category = "SGR1935" if row["source"] == "SGR 1935+2154" else "other"
        key = (year, category)
        rec = groups[key]
        direct, shifted = int(row["narrow_literal"]), int(row["narrow_plus69"])
        rec["n"] += 1
        rec["wins"] += shifted > direct
        rec["losses"] += shifted < direct
        rec["ties"] += shifted == direct
        rec["count_at_label"] += direct
        rec["count_at_plus69"] += shifted
        rec["obs"].add(row["obsid"])
        rec["sources"].add(row["source"])
    assert sum(v["n"] for v in groups.values()) == 944
    assert (groups[("2024", "other")]["n"],
            groups[("2024", "other")]["wins"],
            len(groups[("2024", "other")]["obs"])) == (14, 14, 6)
    assert (groups[("2020", "other")]["wins"],
            groups[("2020", "other")]["n"]) == (64, 64)
    print("Year, target stratum, pairs, +69 wins, ties, losses, distinct visits")
    for key, rec in sorted(groups.items()):
        print(*key, rec["n"], rec["wins"], rec["ties"], rec["losses"],
              len(rec["obs"]), sep=",")
    return groups


def shao_cross_paper(entries):
    by_key = {(r["obsid"], r["number"]): r for r in entries}
    matched = rows(PAIRS)
    assert len(matched) == 10
    seen = set()
    residuals = []
    print("Shao ID, Chu ID, Chu-label minus Shao UTC [s],",
          "Shao UTC minus PR23 +69.184 local UTC [s]")
    for row in matched:
        id_ = row["chu_number"]
        assert id_ not in seen
        seen.add(id_)
        ours = by_key[("5020560106", id_)]
        assert ours["source"] == "SGR 1935+2154"
        shaot = datetime.fromisoformat(row["shao_utc"])
        labelt = datetime.fromisoformat(ours["start"])
        predicted_met = float(ours["predicted_met_tt_s"])
        # The PR23 prediction is local TT seconds relative to the NICER epoch;
        # in calibrated UTC this is 2014 UTC + MET - LEAPINIT. This is the
        # NOMINAL local event time after applying the previously tested shift.
        predicted_utc = UTC_EPOCH + timedelta(
            seconds=predicted_met + TT_MINUS_UTC_2017_PLUS - LEAPINIT)
        label_gap = (labelt - shaot).total_seconds()
        disagreement = (shaot - predicted_utc).total_seconds()
        assert 9.4 < label_gap < 10.2
        assert 60.9 < disagreement < 61.6
        residuals.append(disagreement)
        print(row["shao_id"], id_, f"{label_gap:.3f}",
              f"{disagreement:.3f}", sep=",")
    assert len(residuals) == 10
    assert abs(median(residuals) - 61.275) < 0.002
    print("N=10; median Shao minus frozen PR23 local UTC =",
          f"{median(residuals):.3f} s; range",
          f"{min(residuals):.3f} to {max(residuals):.3f} s")
    return residuals


def negative_controls(entries):
    # These FRBs are frequently confused with the catalog's April 2020
    # and October 2022 burst storms. Time-only claims are NOT justified.
    for day, radio_utc in [
        ("2020-04-28", "2020-04-28T14:34:24"),
        ("2022-10-14", "2022-10-14T19:21:39"),
    ]:
        on_date = [r for r in entries if r["start"].startswith(day)]
        radio = datetime.fromisoformat(radio_utc)
        near = [r for r in on_date
                if abs((datetime.fromisoformat(r["start"]) - radio)
                       .total_seconds()) < 300]
        assert not near, (day, near[:3])
        print("No Chu catalog label within 5 min of", radio_utc,
              "; this does not prove lack of NICER exposure or a counterpart.")


if __name__ == "__main__":
    data = rows(WINDOWS)
    assert len(data) == 1132
    test_time_references()
    year_strata(data)
    shao_cross_paper(data)
    negative_controls(data)
