"""Disentangle time scale (TT/UTC) from barycentric location (local/SSB).

Replays only frozen PR #23 time predictions and the independently published
Shao et al. 2024 UTC labels. Uses Python standard library and no network.
Neither paper's 'UTC' or 'TDB' label is silently reinterpreted as the other.
The "SSB TDB minus local TT" is an APPROXIMATE physical Rømer+relativistic
delay; calendar TDB-TT has a millisecond correction, not calculated here.
"""
import csv
from datetime import datetime, timedelta
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
EPOCH_UTC = datetime(2014, 1, 1)
TT_AT_2014_EPOCH = 67.184  # MJDREFF * 86400
TT_MINUS_UTC_2022 = 69.184
LEAPINIT_2022 = 2.0


def data(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    predicted = {
        (r["obsid"], r["number"]): r
        for r in data(HERE.parent / "catalog-clock" / "windows.csv")
    }
    matches = data(HERE / "shao_2024_table4_pairs.csv")
    assert len(matches) == 10
    rows = []
    for match in matches:
        burst = predicted[("5020560106", match["chu_number"])]
        assert burst["source"] == "SGR 1935+2154"
        met = float(burst["predicted_met_tt_s"])  # PR23 literal local TIME+TIMEZERO
        corrected_met = met + TT_MINUS_UTC_2022
        local_tt = EPOCH_UTC + timedelta(seconds=TT_AT_2014_EPOCH + corrected_met)
        local_utc = EPOCH_UTC + timedelta(seconds=corrected_met - LEAPINIT_2022)
        bary_tdb = datetime.fromisoformat(burst["start"]) + timedelta(
            seconds=TT_MINUS_UTC_2022)
        shao_utc_label = datetime.fromisoformat(match["shao_utc"])
        ssb_minus_local_tt = (bary_tdb - local_tt).total_seconds()
        shao_minus_local_utc = (shao_utc_label - local_utc).total_seconds()
        calendar_difference = (datetime.fromisoformat(burst["start"])
                               - shao_utc_label).total_seconds()
        # This identity is an algebraic *frame bookkeeping* check; it does
        # NOT prove either publication's timestamp is physically correct.
        assert abs(ssb_minus_local_tt - shao_minus_local_utc
                   - calendar_difference) < 0.00001
        rows.append((match["chu_number"], ssb_minus_local_tt,
                     shao_minus_local_utc, calendar_difference))
    assert len(rows) == 10
    record = next(r for r in rows if r[0] == "443")
    assert abs(record[1] - 71.249572) < 0.0001
    assert abs(record[2] - 61.297572) < 0.0001
    print("Chu ID, corrected Chu SSB TDB minus predicted local TT [s],"
          " Shao UTC-labelled time minus predicted local UTC [s],"
          " literal Chu TDB label minus Shao UTC label [s]")
    for item in rows:
        print(item[0], *(f"{v:.6f}" for v in item[1:]), sep=",")
    print("Frame-mixed comparison (Shao - local UTC) median:",
          f"{median(r[2] for r in rows):.6f} seconds")
    print("Do not interpret frame-mixed numbers as a new clock offset.")


if __name__ == "__main__":
    main()
