"""Compare sampled flood-season maxima with OPIDIN annual maximum levels.

Run next to akka_d.csv, dire_d.csv and opidin_maxima.json. Annual flood years
begin in June, as in OPIDIN's historical flood curves. Each OPIDIN station level
has its own local datum. The origin of post-2018 Akka annual values is unverified. WSE differences across years remove fixed offsets.
"""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent


def select(rows, year, coverage, pass_id):
    return [
        row for row in rows
        if f"{year}-06" <= row["time_str"] < f"{year+1}-06"
        and row["time_str"] != "no_data"
        and int(row["reach_q"]) in (0, 1)
        and row["partial_f"] == "0"
        and float(row["obs_frac_n"]) >= coverage
        and int(row["xovr_cal_q"]) <= 1
        and float(row["wse"]) != -999999999999.0
        and (pass_id == "all" or row["pass_id"] == pass_id)
    ]


def main():
    published = {r["year"]: r for r in json.loads((OUT/"opidin_maxima.json").read_text())["data"]}
    fields = ["site", "collection", "year_flood_start", "coverage_min", "pass_filter",
              "eligible_rows", "eligible_peak_window_rows", "opidin_peak_cm_local",
              "sampled_wse_max_m", "peak_wse_u_m", "time_utc", "fraction_on_peak",
              "peak_cycle", "peak_pass", "peak_granule"]
    with (OUT/"peak_matchup.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for site in ("akka", "dire"):
            with (OUT/f"{site}_d.csv").open() as source:
                rows = list(csv.DictReader(source))
            for year in [2023, 2024, 2025]:
                for coverage in [0.8, 0.75, 0.5]:
                    for pass_id in ["all", "363"]:
                        selection = select(rows, year, coverage, pass_id)
                        peak_window = [
                            r for r in selection
                            if f"{year}-10" <= r["time_str"] < f"{year+1}-02"
                        ]
                        data = {
                            "site":site, "collection":"SWOT_L2_HR_RiverSP_D",
                            "year_flood_start": year, "coverage_min": coverage,
                            "pass_filter": pass_id, "eligible_rows": len(selection),
                            "eligible_peak_window_rows":len(peak_window),
                            "opidin_peak_cm_local": published[year][f"level_{site}"],
                        }
                        # Restrict reporting to years with data during plausible peak period.
                        if peak_window:
                            peak = max(selection, key=lambda r: float(r["wse"]))
                            data.update(sampled_wse_max_m=peak["wse"],
                                        peak_wse_u_m=peak["wse_u"],
                                        time_utc=peak["time_str"], fraction_on_peak=peak["obs_frac_n"],
                                        peak_cycle=peak["cycle_id"], peak_pass=peak["pass_id"],
                                        peak_granule=peak["granuleUR"])
                        writer.writerow(data)


if __name__ == "__main__":
    main()
