"""Rebuild the optical outburst chronology from the frozen epoch rows."""
from __future__ import annotations

import csv
from collections import defaultdict
from datetime import datetime, timedelta, timezone
import json
import math
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent
MJD_ZERO = datetime(1858, 11, 17, tzinfo=timezone.utc)


def read(name):
    with (ROOT / "data" / name).open() as stream:
        return list(csv.DictReader(stream))


def date(mjd):
    return (MJD_ZERO + timedelta(days=mjd)).date().isoformat()


def separation(a, b):
    dx = (a[0] - b[0]) * math.cos(math.radians((a[1] + b[1]) / 2))
    return 3600 * math.hypot(dx, a[1] - b[1])


def centroid(rows, ra, dec):
    return [median(float(r[c]) for r in rows) for c in (ra, dec)]


def make_episodes(observations, gap_days):
    groups = []
    for item in sorted(observations):
        if not groups or item[0] - groups[-1][-1][0] > gap_days:
            groups.append([])
        groups[-1].append(item)
    return groups


def main():
    ptf = [r for r in read("ptf.csv") if r["goodflag"] == "1" and r["qaflags"] == "0"]
    direct = [r for r in read("ztf_direct.csv") if int(r["catflags"]) == 0 and float(r["magerr"]) <= 0.25]
    alerts = [r for r in read("ztf_alerts.csv") if int(r["isdiffpos"]) == 1 and float(r["sigmapsf"]) <= 0.25 and float(r["rb"]) >= 0.55]
    points = [(float(r["obsmjd"]), float(r["mag_autocorr"]), "PTF", {"1":"g","2":"R"}[r["fid"]]) for r in ptf]
    points += [(float(r["mjd"]), float(r["mag"]), "ZTF direct", r["filtercode"][1]) for r in direct]
    points += [(float(r["mjd"]), float(r["magpsf"]), "ZTF difference", {"1":"g","2":"r","3":"i"}[r["fid"]]) for r in alerts]
    bright = [p for p in points if p[1] <= 20.0]
    limits = [r for r in read("ztf_limits.csv") if float(r["diffmaglim"]) >= 20.3]
    episodes = []
    for group in make_episodes(bright, 45):
        start, end = group[0][0], group[-1][0]
        nights = len({int(p[0]) for p in group})
        pre = [float(r["mjd"]) for r in limits if start-70 <= float(r["mjd"]) < start]
        by_survey = defaultdict(list)
        for _, mag, survey, band in group:
            by_survey[(survey, band)].append(mag)
        episodes.append({
            "first_utc": date(start), "last_bright_utc": date(end),
            "first_mjd": round(start, 5), "last_bright_mjd": round(end, 5),
            "bright_utc_nights": nights, "instruments": sorted({p[2] for p in group}),
            "brightest_sample_by_product_band": {s+" "+b: round(min(m), 3) for (s, b), m in sorted(by_survey.items())},
            "previous_difference_limit_within_70d": date(max(pre)) if pre else None})
    coord = centroid(alerts, "ra", "dec")
    position = {}
    for year in (2011, 2014):
        group = [r for r in ptf if date(float(r["obsmjd"])).startswith(str(year))]
        p = centroid(group, "s_ra", "s_dec")
        position[str(year)] = {"n": len(group), "median_radec_deg": p, "separation_from_ztf_arcsec": round(separation(coord, p), 3)}
    nights = defaultdict(list)
    for r in direct:
        if r["filtercode"] == "zg":
            nights[int(float(r["mjd"]))].append((float(r["mjd"]), float(r["mag"]), float(r["magerr"])))
    for r in ptf:
        if r["fid"] == "1":
            nights[int(float(r["obsmjd"]))].append((float(r["obsmjd"]), float(r["mag_autocorr"]), float(r["magerr_auto"])))
    night_g = {n: [median(col) for col in zip(*a)] for n, a in sorted(nights.items())}
    spans = {}
    for name, start, end in [("2014 PTF g",56949,56955),("2018 ZTF g",58426,58434),("2020 ZTF g",59182,59204),("2023 ZTF g",60141,60148),("2024 ZTF g",60511,60516)]:
        s,e = night_g[start],night_g[end]
        spans[name] = {"start_median_mag": round(s[1],3), "end_median_mag": round(e[1],3),
                       "sampled_span_days": round(e[0]-s[0], 3), "change_mag": round(e[1]-s[1],3)}
    # Under the same bright cut, contrast recurrent high states in different seasons.
    # The comparison uses consecutive intervals of sampled bright points, without
    # interpolating across gaps or assuming a measured true peak.
    peak_bounds = {}
    for y in [2018,2020,2023,2024]:
        g_nights = [v for v in night_g.values() if date(v[0]).startswith(str(y)) and v[1] <= 19.1]
        if g_nights:
            peak_bounds[str(y)]=round(max(v[0] for v in g_nights)-min(v[0] for v in g_nights),3)
    field = read("ps1_neighbours.csv")
    near = sorted((separation(coord,(float(x["raMean"]),float(x["decMean"]))),x["objID"]) for x in field)
    target2014=[r for r in ptf if date(float(r["obsmjd"])).startswith("2014")]
    out = {"ptf_good_epochs":len(ptf),"ptf_2014_nights":len({int(float(r["obsmjd"])) for r in target2014}),
           "ztf_direct_good_epochs":len(direct),"ztf_difference_good_epochs":len(alerts),
           "ztf_median_radec_deg":coord,"ptf_positions":position,"bright_episodes":episodes,
           "number_of_multinight_bright_episodes":sum(e["bright_utc_nights"]>1 for e in episodes),
           "episode_counts_at_gap_days":{str(t):len(make_episodes(bright,t)) for t in (25,45,70,90)},
           "sampled_g_declines":spans,"sampled_g_le_19p1_separations_days":peak_bounds,
           "ps1_two_or_more_detections_within_120_arcsec":len(field),
           "nearest_ps1_objects":[{"objID":i,"arcsec":round(s,3)} for s,i in near[:2]]}
    (ROOT/"measurements.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
