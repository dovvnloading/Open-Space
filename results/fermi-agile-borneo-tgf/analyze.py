"""Reproduce the AGILE/Fermi catalog association, geometry, and GBM TTE check."""
import argparse
import csv
import hashlib
import json
import re
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
from astropy.io import fits

GBM = "https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/tgf/"
TTE = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/2016/03/20/current/"
URLS = {
    "off.csv": GBM + "gbm_tgf_catalog_offline.csv",
    "ww.csv": GBM + "gbm_tgf_catalog_wwlln.csv",
    "agile.html": "https://www.ssdc.asi.it/mcal3tgfcat/",
    "radio_event.txt": ("https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/tgf/cat_files/"
                       "tgf_wwlln_data/gbm_WWLLN.oTGF160320330.txt"),
    **{f"{d}.gz": TTE + f"glg_tte_{d}_160320_07z_v00.fit.gz" for d in ("b0", "b1")},
}
FERMI_ID = "oTGF160320330"
EPOCH = datetime(2000, 1, 1)


def utc_microseconds(t):
    delta = t - EPOCH
    return (delta.days * 86400 + delta.seconds) * 1_000_000 + delta.microseconds


def read_sources(cache):
    cache.mkdir(parents=True, exist_ok=True)
    sources = {}
    for name, url in URLS.items():
        path = cache / name
        if not path.exists():
            tmp = path.with_suffix(".partial")
            with urllib.request.urlopen(url, timeout=120) as response:
                tmp.write_bytes(response.read())
            tmp.replace(path)
        sources[name] = {"url": url, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return sources


def crossmatch(cache):
    fermi = list(csv.DictReader((cache / "off.csv").read_text().splitlines(), skipinitialspace=True))
    radio = list(csv.DictReader((cache / "ww.csv").read_text().splitlines(), skipinitialspace=True))
    # The full-resolution UTC is embedded in the AGILE new ID (the display UTC can be truncated).
    # Part II selects gamma-ray data without imposing a WWLLN coincidence. Part I rows can
    # duplicate these same photons and are used exclusively for the radio check below.
    agile = [next(csv.reader([r])) for r in re.findall(
        r'this\[\d+\]\s*=\s*new catEntry\s*\(([^\n]*?)\);', (cache / "agile.html").read_text())]
    photon = {}
    for row in agile:
        if row[22] in {"5", "7"}:
            t = utc_microseconds(datetime.strptime(row[18], "%Y%m%d.%H%M%S.%f"))
            photon.setdefault(t, row)
    gbm = sorted((utc_microseconds(datetime.fromisoformat(r["Date"] + "T" + r["UTC"])), r) for r in fermi)
    at = np.array(sorted(photon))
    ft = np.array([x[0] for x in gbm])
    lo, hi = max(at[0], ft[0]), min(at[-1], ft[-1])
    idx = np.where((ft >= lo) & (ft <= hi))[0]
    ia = (at >= lo) & (at <= hi)
    common = at[ia]

    def nearest(times):
        pos = np.searchsorted(common, times)
        left = np.maximum(0, pos - 1)
        right = np.minimum(len(common) - 1, pos)
        near = np.where(abs(common[left] - times) < abs(common[right] - times), left, right)
        return common[near], abs(common[near] - times)

    match, distance = nearest(ft[idx])
    ranked = np.argsort(distance)
    top = []
    for row in ranked[:5]:
        f = gbm[idx[row]][1]
        a = photon[match[row]]
        top.append({"fermi_id": f["#OS_ID"], "agile_utc_id": a[18], "gap_s": distance[row].item() / 1e6,
                    "nadir_distance_km": round(surface(float(f["Lat"]), float(f["Lon"]), float(a[3]), float(a[2])), 2)})
    best = ranked[0]
    fer = gbm[idx[best]][1]
    ag = photon[match[best]]
    assert fer["#OS_ID"] == FERMI_ID, "The current catalogs no longer select the published pair."
    ar = [r for r in agile if r[6] == ag[6] and r[22] == "3" and r[16]]
    fr = [r for r in radio if r["#ID"].strip() == FERMI_ID]
    assert len(ar) == len(fr) == 1
    # Common shift for all AGILE times (equivalently shift GBM in the opposite
    # direction). Recompute match count for the whole search on each shift. Keep
    # edge margins fixed and 1,000 nonzero offsets in each direction.
    margins = idx[(ft[idx] >= lo + 1_000_100_000) & (ft[idx] <= hi - 1_000_100_000)]
    shifts = np.r_[np.arange(-1000, 0), np.arange(1, 1001)]
    controls = np.array([np.sum(nearest(ft[margins] + shift * 1_000_000)[1] <= 100_000) for shift in shifts])
    observed = int(np.sum(nearest(ft[margins])[1] <= 100_000))
    return fer, ag, fr[0], ar[0], {
        "fermi_full": len(fermi), "agile_photon_rows": sum(r[22] in {"5", "7"} for r in agile),
        "agile_unique_photon_times": len(photon), "fermi_common_interval": len(idx),
        "agile_unique_common_interval": int(ia.sum()), "nearest_pairs": top,
        "fermi_events_with_match_within_100ms": int(sum(distance <= 100_000)), "time_shift_window_s": 0.1,
        "shift_integers_s": [-1000, -1, 1, 1000], "edge_matched_fermi": len(margins),
        "observed_edge_matched": observed, "shifts_with_at_least_one_pair": int(sum(controls > 0)),
        "total_pairs_in_2000_shifts": int(controls.sum()),
    }


def xyz(lat, lon, h):
    lat, lon = np.deg2rad([lat, lon])
    f = 1 / 298.257223563
    e2 = f * (2 - f)
    n = 6378.137 / np.sqrt(1 - e2 * np.sin(lat)**2)
    return np.array([(n + h)*np.cos(lat)*np.cos(lon), (n + h)*np.cos(lat)*np.sin(lon),
                     (n * (1 - e2) + h)*np.sin(lat)])


def surface(lat0, lon0, lat1, lon1):
    p, q, dl = np.deg2rad([lat0, lat1, lon1 - lon0])
    return 6371.0 * 2 * np.arcsin(np.sqrt(np.sin((q-p)/2)**2 + np.cos(q)*np.cos(p)*np.sin(dl/2)**2))


def angle(u, v):
    return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(u, v)), np.dot(u, v)))


def geometry(f, a, fr, ar):
    lat, lon = float(fr["WWLLN_Lat"]), float(fr["WWLLN_Lon"])
    sat_f = xyz(float(f["Lat"]), float(f["Lon"]), float(f["Alt"]))
    sat_a = xyz(float(a[3]), float(a[2]), float(a[10]))
    src = xyz(lat, lon, 15)
    u, v = sat_f - src, sat_a - src
    zenith = xyz(lat, lon, 16) - src
    # Geodesically displace source by up to 15 km (spherical direct formula).
    angles = []
    p, l = np.deg2rad([lat, lon])
    for h in range(10, 21):
        for km in (0, 5, 10, 15):
            d = km / 6371
            for azi in np.deg2rad(np.arange(0, 360, 5)):
                q = np.arcsin(np.sin(p)*np.cos(d) + np.cos(p)*np.sin(d)*np.cos(azi))
                z = l + np.arctan2(np.sin(azi)*np.sin(d)*np.cos(p), np.cos(d) - np.sin(p)*np.sin(q))
                s = xyz(np.rad2deg(q), np.rad2deg(z), h)
                angles.append(angle(sat_f-s, sat_a-s))
    delay_ms = (np.linalg.norm(u) - np.linalg.norm(v)) / 299.792458
    return {
        "wwlln_lat_lon_deg": [lat, lon], "agile_radio_lat_lon_deg": [float(ar[15]), float(ar[14])],
        "assumed_altitude_km": 15, "fermi_xyz_wgs84_km": list(sat_f), "agile_xyz_wgs84_km": list(sat_a),
        "spacecraft_distance_km": np.linalg.norm(sat_f-sat_a),
        "fermi_agile_slant_km": [np.linalg.norm(u), np.linalg.norm(v)],
        "fermi_agile_off_zenith_deg": [angle(u, zenith), angle(v, zenith)],
        "separation_deg": angle(u, v), "position_height_grid_separation_deg": [min(angles), max(angles)],
        "fermi_minus_agile_propagation_ms": delay_ms,
        "fermi_wwlln_offset_ms_catalog": float(fr["Separation_ms"]),
        "agile_wwlln_offset_ms_catalog": float(ar[16])/1000,
    }


def gbm_photons(cache, f):
    photons, background = [], {}
    width = float(f["Width_ms"]) / 1000
    for det in ("b0", "b1"):
        with fits.open(cache / f"{det}.gz", memmap=False) as hdus:
            eb = hdus["EBOUNDS"].data
            ev = hdus["EVENTS"].data
            t = ev["TIME"] - float(f["MET"])
            channel = ev["PHA"]
            # Equal-duration sidebands on each side, with 10 ms guard around the pulse.
            side = (abs(t) >= 0.01) & (abs(t) < 1.01)
            background[det] = {"sideband_counts_2s_all": int(side.sum()),
                               "expected_all_in_window": float(side.sum()*width/2)}
            for i in np.flatnonzero((t >= 0) & (t < width)):
                ch = int(channel[i])
                photons.append({"detector": det, "dt_us": round(t[i]*1e6, 2), "pha": ch,
                                "emin_kev": round(float(eb["E_MIN"][ch]), 2), "emax_kev": round(float(eb["E_MAX"][ch]), 2)})
    return sorted(photons, key=lambda r: r["dt_us"]), background


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", required=True, type=Path, help="Input downloads, about 60 MB total")
    p.add_argument("--output", type=Path, default=Path(__file__).parent)
    args = p.parse_args()
    sources = read_sources(args.cache)
    f, a, fr, ar, scan = crossmatch(args.cache)
    geom = geometry(f, a, fr, ar)
    photons, background = gbm_photons(args.cache, f)
    event = {"fermi": f, "fermi_wwlln": fr, "agile_selection": {
        "utc_id": a[18], "lon_deg": float(a[2]), "lat_deg": float(a[3]), "alt_km": float(a[10]),
        "orbit": int(a[4]), "ml_counts": float(a[12]), "ml_error": float(a[13]),
        "t50_ms": float(a[7]), "t50_error_ms": float(a[8]), "event_type": int(a[19]),
        "part_ii_flag": int(a[22]), "part_i_wwlln_drift_flag": int(ar[22])}}
    lines = (args.cache / "radio_event.txt").read_text().splitlines()
    associations, n_strokes = map(int, lines[4].split())
    strokes = sorted((abs(float(row.split()[3])), row.split()) for row in lines[5 + associations:])
    assert len(strokes) == n_strokes and float(strokes[0][1][1]) == float(fr["WWLLN_Lon"])
    field = {"n_strokes_in_20_minutes": n_strokes, "closest_abs_offset_s": strokes[0][0],
             "next_closest_abs_offset_s": strokes[1][0]}
    report = {"event": event, "scan": scan, "geometry": geom, "radio_field": field,
              "tte_background": background, "sources": sources}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "evidence.json").write_text(json.dumps(report, indent=2) + "\n")
    with (args.output / "gbm_photons.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=photons[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(photons)


if __name__ == "__main__":
    main()
