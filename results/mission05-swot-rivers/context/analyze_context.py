"""Summarize the published gauge and CHIRPS extracts without gap-filling."""
import calendar
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean

HERE = Path(__file__).resolve().parent


def rows(name):
    with (HERE / name).open() as source:
        return list(csv.DictReader(source))


def main():
    rain = defaultdict(dict)
    for row in rows("chirps_wet_season_dekads.csv"):
        key = row["PCODE"], row["date"][:7]
        assert row["version"] == "final" and row["date"] not in rain[key]
        rain[key][row["date"]] = float(row["rfh"])
    monthly = {}
    for key, observations in rain.items():
        assert sorted(day[-2:] for day in observations) == ["01", "11", "21"]
        monthly[key] = sum(observations.values())

    for code in sorted({key[0] for key in monthly}):
        february = {y: monthly[code, f"{y}-02"] for y in range(1982, 2026)}
        normal = mean(february[y] for y in range(1991, 2021))
        # Annual wet-season index is the year containing January.
        season = {}
        for year in range(1982, 2026):
            spans = [(year-1, 10), (year-1, 11), (year-1, 12), (year, 1), (year, 2), (year, 3)]
            season[year] = sum(monthly[code, f"{yr}-{mo:02}"] for yr, mo in spans)
        seasonal_mean = mean(season[y] for y in range(1991, 2021))
        print(code, "February precipitation mm, climatology/2023/2024/2025",
              *(round(x, 1) for x in (normal, february[2023], february[2024], february[2025])))
        feb_rank = 1 + sum(monthly[code, f"{y}-02"] < february[2024] for y in range(1981, 2026))
        print("  February 2024 rank among 45 years", feb_rank)
        print("  Oct–Mar precipitation mm, climatology/2023/2024/2025",
              *(round(x, 1) for x in (seasonal_mean, season[2023], season[2024], season[2025])))

    flows = defaultdict(dict)
    for row in rows("zra_daily_gauges.csv"):
        flows[row["station"]][row["date"]] = float(row["q_m3s"])
    normal = {int(r["calendar_month"]): float(r["mean_q_m3s"])
              for r in rows("chavuma_historical_seasonal_flow.csv")}
    for month in (2, 3, 4, 5):
        key = f"2024-{month:02}-"
        values = [v for day, v in flows["Chavuma"].items() if day.startswith(key)]
        assert len(values) == calendar.monthrange(2024, month)[1]
        print("Chavuma", key[:7], "gauge mean, historic gauge mean, difference %",
              round(mean(values), 1), normal[month], round(100*(mean(values)/normal[month] - 1), 1))

    for station, obs in flows.items():
        for year in (2023, 2024):
            keys = [f"04-{d:02}" for d in range(1, 31)
                    if f"{year}-04-{d:02}" in obs and f"{year+1}-04-{d:02}" in obs]
            x, y = [mean(obs[f"{yr}-{d}"] for d in keys) for yr in (year,year+1)]
            print(station, "April paired", year, year+1, len(keys), round(x,1), round(y,1),
                  round(100*(y/x-1),1))


if __name__ == "__main__":
    main()
