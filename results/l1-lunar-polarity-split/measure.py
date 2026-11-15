"""Measure the held-out lunar interval and L1 controls from the one-minute evidence."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
d = pd.concat([pd.read_csv(ROOT / name, parse_dates=["utc"], index_col="utc")
               for name in ["minute.csv.gz", "orbits.csv.gz"]], axis=1)


def interval(start, end):
    """Select UTC minute starts in a left-closed, right-open interval on 2019-01-31."""
    lo = pd.Timestamp(f"2019-01-31 {start}", tz="UTC")
    hi = pd.Timestamp(f"2019-01-31 {end}", tz="UTC")
    return d.loc[(d.index >= lo) & (d.index < hi)]


def magnetic(window, probes):
    result = {}
    for probe in probes:
        values = window[probe + "_bz"].dropna()
        result[probe] = {"minutes": len(values), "mean_nT": round(values.mean(), 3),
                         "northward_minutes": int((values > 0).sum())}
    return result


up = interval("17:30", "18:30")
# The hour after the selected L1 hour encloses the ~42–48 min radial travel delays.
down = interval("18:30", "19:30")
re = 6378
pos = {k: up[[k+"_x_km", k+"_y_km", k+"_z_km"]].mean().to_numpy() for k in ["ace", "wind", "dscovr"]}
pos.update({k: down[[k+"_x_km", k+"_y_km", k+"_z_km"]].mean().to_numpy() for k in ["thb", "thc"]})
vx = float(abs(up.ace_vx.median()))
transit = {k: round((pos[k][0] - pos["thb"][0]) / vx / 60, 2) for k in ["ace", "wind", "dscovr"]}
lag_hours = {}
for k in ["ace", "wind", "dscovr"]:
    means = [d[k+"_bz"].shift(lag).loc[up.index].mean() for lag in range(-30, 31)]
    lag_hours[k] = [round(min(means), 3), round(max(means), 3)]
# Compare independently on valid, common minutes; default study window and its shoulders.
joined = down[["thb_bz", "thc_bz", "omni_bz"]].dropna()
sources = down.loc[joined.index, "omni_source"].value_counts().to_dict()
omni_shift = {}
for k in range(-15, 16):
    shifted = d.omni_bz.shift(k).loc[down.index]
    if shifted.count() >= 40:
        omni_shift[k] = round(shifted.mean(), 3)
result = {
    "upstream_1730_1830": magnetic(up, ["ace", "wind", "dscovr"]),
    "lunar_1830_1930": magnetic(down, ["thb", "thc", "omni"]),
    "lunar_both_north_omni_south_common_minutes": int(((joined.thb_bz>0)&(joined.thc_bz>0)&(joined.omni_bz<0)).sum()),
    "lunar_omni_common_minutes": len(joined), "omni_imf_source_ids_in_comparison": sources,
    "median_abs_ace_vx_km_s": round(vx, 3), "nominal_radial_minutes_to_thb": transit,
    "ace_wind_radial_delay_minutes": round((pos["wind"][0]-pos["ace"][0])/vx/60, 2),
    "ace_wind_transverse_km": round(np.linalg.norm((pos["wind"]-pos["ace"])[1:])),
    "wind_dscovr_transverse_km": round(np.linalg.norm((pos["wind"]-pos["dscovr"])[1:])),
    "lunar_probe_separation_RE": round(np.linalg.norm(pos["thb"]-pos["thc"])/re,3),
    "mean_positions_GSE_RE": {k:np.round(v/re, 3).tolist() for k,v in pos.items()},
    "hour_mean_bz_shift_sensitivity_min_max_nT_at_L1_30min":lag_hours,
    "omni_hour_mean_with_15min_shifts_min_max_nT":[min(omni_shift.values()), max(omni_shift.values())],
    "ace_rotation_check_rms_nT":round(float(np.sqrt(np.mean(d.ace_rotation_residual.dropna()**2))),5),
    "lunar_plasma_good_bins":{k:int(down[k+"_vx"].count()) for k in ["thb","thc"]},
    "lunar_median_vx_km_s":{k:round(down[k+"_vx"].median(),2) for k in ["thb","thc"]},
    "common_later_southward_control_2000_2030":magnetic(interval("20:00","20:30"),["ace","wind","dscovr"]),
    "lunar_1815_1915_alternate_window": magnetic(interval("18:15", "19:15"), ["thb","thc","omni"]),
}
print(json.dumps(result, indent=2))
