"""Recover AAOT optical, seasonal, hydrological and weather evidence for November 2025."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.parse import urlencode

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
AOC = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?"
BASE = "https://marine-api.open-meteo.com/v1/marine?latitude=45.3139&longitude=12.5083&start_date=2025-11-01&end_date=2025-11-10&"
AIR_URL = ("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=LIPZ&data=all&year1=2025&month1=11"
           "&day1=1&year2=2025&month2=11&day2=10&tz=Etc%2FUTC&format=onlycomma&latlon=yes&elev=no"
           "&missing=M&trace=T&direct=no&report_type=3&report_type=4")
PO = "https://www.arpae.it/it/temi-ambientali/meteo/report-meteo/bollettini-mensili/bm_202511.pdf"
PDF_DIGEST = "742520a217efab4ef234cdc39fff1aaeb3cab02c77e103d52a84525e33605e70"
OCI = "https://pace.oceansciences.org/data/matchup_Rrs.csv"


def dataset(cache, manifest, name, url):
    path = cache / name
    if not path.exists():
        data = subprocess.check_output(["curl", "-fsSL", "--max-time", "120", "--retry", "2", url])
        if data[:50].lower().strip().startswith((b"<html", b"<!doctype")):
            raise ValueError(f"Expected data, received HTML for {url}")
        path.write_bytes(data)
    manifest.append({"file": name, "url": url, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    return path


def daily_curve(frame):
    channels = frame[[f"Lwn_f/Q[{w}nm]" for w in (560, 620, 667)]]
    x = np.log(channels.where(channels > 0))
    return x.iloc[:, 1] - (47 * x.iloc[:, 0] + 60 * x.iloc[:, 2]) / 107


def main(cache, out, matchup):
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    if matchup is None:
        matchup = dataset(cache, manifest, "matchup_Rrs.csv", OCI)
    else:
        manifest.append({"file": matchup.name, "url": OCI, "sha256": hashlib.sha256(matchup.read_bytes()).hexdigest()})
    data = pd.read_csv(matchup, comment="#", na_values=-999, low_memory=False)
    pair = data.loc[(data.FIELD_cruise == "aoc_v4l15_AAOT") &
                    data.SAT_datetime.str.startswith(("2025-11-06", "2025-11-07"))].sort_values("SAT_datetime")
    if len(pair) != 2:
        raise ValueError("Expected one released AAOT matchup on each day")
    keys = ["FIELD_datetime", "SAT_datetime", "SAT_parent_file", "FIELD_latitude", "FIELD_longitude",
            "SAT_cv", "SAT_valid_pixels", "SAT_QWIP_Score", "SAT_senz", "SAT_solz", "SAT_elev", "SATFIELD_time_diff"]
    pair[keys].to_csv(out / "matchup_quality.csv", index=False)
    bands = []
    for _, obs in pair.iterrows():
        for w in sorted(set(range(400, 701, 5)) | {412, 443, 532}):
            row = {"date": obs.SAT_datetime[:10], "wavelength_nm": w}
            for src, dst in [("SAT_Rrs", "oci_rrs"), ("FIELD_Rrs", "field_rrs"),
                             ("SAT_Rrs_unc", "oci_unc"), ("FIELD_Rrs_unc", "field_unc"),
                             ("SAT_Rrs_spatstd", "oci_spatial_sd")]:
                row[dst] = obs[f"{src}_{w}"]
            bands.append(row)
    pd.DataFrame(bands).to_csv(out / "spectra.csv", index=False)
    stack = []
    for year in range(2017, 2026):
        query = urlencode(dict(site="AAOT", year=year, month=10, day=7, year2=year, month2=12,
                               day2=7, LWN15=1, AVG=20, if_no_html=1))
        path = dataset(cache, manifest, f"aaot_l15_{year}.csv", AOC + query)
        raw = pd.read_csv(path, skiprows=5, na_values=-999)
        stack.append(pd.DataFrame({"date": pd.to_datetime(raw["Date(dd-mm-yyyy)"], format="%d:%m:%Y"),
                                   "C": daily_curve(raw), "Lwn560": raw["Lwn_f/Q[560nm]"]}))
    daily = pd.concat(stack).set_index("date").sort_index()
    days = daily.C
    brightness_rows = []
    past_days = daily.loc[daily.index.year < 2025]
    for endpoint in ("2025-11-06", "2025-11-07"):
        current = daily.loc[endpoint]
        for tolerance in (0.1, 0.15, 0.2):
            comparable = past_days.loc[(past_days.Lwn560 / current.Lwn560 - 1).abs() <= tolerance, "C"].dropna()
            brightness_rows.append({"endpoint": endpoint, "fractional_tolerance": tolerance,
                                    "daily_C": current.C, "n_prior_days": len(comparable),
                                    "historical_C_median": comparable.median(),
                                    "historical_C_lower_or_equal": int((comparable <= current.C).sum())})
    pd.DataFrame(brightness_rows).to_csv(out / "green_brightness_controls.csv", index=False)
    increments = (days - days.shift(freq="1D")).dropna()
    increments.to_csv(out / "field_daily_increments.csv", header=["delta_C560_620_667"])
    match = increments.loc["2025-11-07"]
    past = increments.loc[increments.index.year < 2025]
    print(f"Observed {match:.6f}; {sum(past <= match)} of {len(past)} previous matched-season days <= event.")
    days.loc["2025-11-03":"2025-11-12"].to_csv(out / "field_daily_event.csv", header=["C560_620_667"])
    wave_url = BASE + "models=meteofrance_wave&hourly=wave_height,wave_period,wave_direction&timezone=UTC"
    current_url = BASE + "hourly=ocean_current_velocity,ocean_current_direction,sea_level_height_msl&timezone=UTC"
    wave = json.loads(dataset(cache, manifest, "wave_explicit.json", wave_url).read_text())
    current = json.loads(dataset(cache, manifest, "current.json", current_url).read_text())
    frame = pd.merge(pd.DataFrame(wave["hourly"]), pd.DataFrame(current["hourly"]), on="time")
    frame["wave_latitude"] = wave["latitude"]
    frame["wave_longitude"] = wave["longitude"]
    frame["current_latitude"] = current["latitude"]
    frame["current_longitude"] = current["longitude"]
    frame.to_csv(out / "modeled_marine_hours.csv", index=False)
    air = pd.read_csv(dataset(cache, manifest, "wind.csv", AIR_URL), na_values="M")
    cols = ["station", "valid", "lon", "lat", "drct", "sknt", "gust", "metar"]
    air.loc[air.valid.str.startswith(("2025-11-06", "2025-11-07")), cols].to_csv(out / "airport_wind.csv", index=False)
    pdf = dataset(cache, manifest, "po.pdf", PO)
    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_DIGEST:
        raise ValueError("Po bulletin changed. Recheck the printed values on p. 37.")
    # Visual transcription of the Pontelagoscuro column of Table 1. The PDF embeds the table as an image.
    pd.DataFrame({"date": pd.date_range("2025-11-01", periods=12).date,
                  "pontelagoscuro_daily_m3_s": [1200,1201,1136,1304,1422,1241,1118,1068,1035,1006,978,959]}
                 ).to_csv(out / "po_downstream_gauge.csv", index=False)
    pd.DataFrame(manifest).to_csv(out / "sources.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--matchup", type=Path)
    parser.add_argument("--out", type=Path, default=ROOT)
    args = parser.parse_args()
    main(args.cache, args.out, args.matchup)
