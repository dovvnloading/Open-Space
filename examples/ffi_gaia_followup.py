"""Test resolved Gaia DR3 photometry around recurrent TESS FFI centroids."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from astropy.timeseries import LombScargle

from open_space.tess_archive import atomic_bytes, download
from open_space.tess_eb_catalog import separation

TAP = "https://gea.esac.esa.int/tap-server/tap/sync?"
DATALINK = "https://gea.esac.esa.int/data-server/data?"
RADIUS_ARCSEC = 30
PERIOD_RANGE_DAYS = (0.15, 8.0)


def receipt(url, cache):
    path = cache / (hashlib.sha256(url.encode()).hexdigest() + ".csv")
    if not path.exists():
        atomic_bytes(path, download(url, 8_000_000))
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if "source_id" not in (reader.fieldnames or ()):
        raise ValueError("Unrecognized Gaia CSV response")
    return list(reader), {"url": url, "sha256": hashlib.sha256(raw).hexdigest()}


def epoch_period(gaia_id, cache):
    url = DATALINK + urlencode({"RETRIEVAL_TYPE": "EPOCH_PHOTOMETRY", "ID": gaia_id,
                                "FORMAT": "CSV", "BAND": "G", "RELEASE": "Gaia DR3"})
    rows, provenance = receipt(url, cache)
    if any(r["source_id"] != gaia_id for r in rows):
        raise ValueError("Gaia epoch series did not match requested source")
    usable = [r for r in rows if r["variability_flag_g_reject"] == "false"
              and r["rejected_by_photometry"] == "false" and r["g_transit_flux"]
              and r["g_transit_flux_error"]]
    t, f, err = [np.asarray([float(r[k]) for r in usable]) for k in
                 ("g_transit_time", "g_transit_flux", "g_transit_flux_error")]
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(err) & (err > 0)
    t, f, err = t[ok], f[ok], err[ok]
    if len(t) < 20 or np.ptp(t) < 100:
        return {"status": "insufficient_epochs", "usable_epochs": len(t), "receipt": provenance}
    t = t - np.min(t)
    ls = LombScargle(t, f, err)
    freq, power = ls.autopower(minimum_frequency=1/PERIOD_RANGE_DAYS[1],
                               maximum_frequency=1/PERIOD_RANGE_DAYS[0], samples_per_peak=12)
    peak = int(np.argmax(power))
    model = ls.model(t, freq[peak])
    p0 = np.average(f, weights=1/err**2)
    improvement = 1 - np.sum(((f-model)/err)**2) / np.sum(((f-p0)/err)**2)
    if not 0 <= improvement <= 1:
        raise ValueError("Unexpected weighted variance fraction")
    return {"status": "measured", "usable_epochs": len(t), "span_days": float(np.ptp(t)),
            "blind_period_days": float(1 / freq[peak]),
            "weighted_variance_fraction": float(improvement), "frequencies_tested": len(freq),
            "receipt": provenance}


def case(source, field, cache):
    ra, dec = source["ra_deg"], source["dec_deg"]
    select = ("SELECT source_id,ra,dec,phot_g_mean_mag,phot_variable_flag,has_epoch_photometry "
              "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
              f"CIRCLE('ICRS',{ra:.9f},{dec:.9f},{RADIUS_ARCSEC/3600:.10f})) ORDER BY source_id")
    rows, provenance = receipt(TAP + urlencode({"REQUEST": "doQuery", "LANG":"ADQL",
                                                "FORMAT":"csv", "QUERY": select}), cache)
    if len(rows) >= 2000:
        raise ValueError("Sky-cone default row limit may have been reached")
    p1, p2 = (source[k]["period_days"] for k in ("discovery", "validation"))
    baselines = [pr["stop_btjd"] - pr["start_btjd"] for pr in field["products"]]
    neighbors = []
    for row in rows:
        if row["has_epoch_photometry"].lower() != "true":
            continue
        ep = epoch_period(row["source_id"], cache)
        distance = separation(ra, dec, float(row["ra"]), float(row["dec"]))
        comparison = None
        if ep["status"] == "measured":
            freq = 1/ep["blind_period_days"]
            # Frequency offsets in units of each TESS sector's Fourier resolution.
            deltas = [abs(freq - 1/p) * span for p,span in zip((p1,p2),baselines)]
            comparison = {"frequency_offsets_rayleigh": deltas,
                          "meets_direct_match": (distance <= 15 and max(deltas) <= 1
                                                 and ep["weighted_variance_fraction"] >= .4)}
        neighbors.append({"source_id": row["source_id"], "ra_deg": float(row["ra"]),
                          "dec_deg": float(row["dec"]), "g_mag": float(row["phot_g_mean_mag"]),
                          "separation_arcsec": distance, "resolved_epoch_search": ep,
                          "comparison": comparison})
    neighbors.sort(key=lambda n: n["separation_arcsec"])
    return {"field_id": field["field_id"], "source_id": source["source_id"],
            "centroid_deg": [ra, dec], "tess_periods_days": [p1,p2],
            "sector_baselines_days": baselines, "cataloged_stars_in_cone": len(rows),
            "count_with_epochs": len(neighbors), "count_marked_variable": sum(
                r["phot_variable_flag"] == "VARIABLE" for r in rows),
            "gaia_epoch_series": neighbors, "gaia_cone_receipt": provenance}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field_result", type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if not 1 <= args.workers <= 10:
        parser.error("Choose 1 to 10 Gaia query workers")
    raw = args.field_result.read_bytes()
    fields = json.loads(raw)["fields"]
    targets = [(s, field) for field in fields for s in field["sources"] if s["candidate"]]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(lambda sf: case(*sf, args.cache), targets))
    print(json.dumps({"ffi_result_sha256": hashlib.sha256(raw).hexdigest(),
                      "archive": "ESA Gaia DR3", "cone_radius_arcsec": RADIUS_ARCSEC,
                      "gaia_period_search_days": list(PERIOD_RANGE_DAYS),
                      "method": "Weighted floating-mean single-sinusoid, G-band epoch photometry",
                      "direct_match_thresholds": {"separation_arcsec_max": 15,
                                                  "fraction_of_rayleigh_max_both_tess_sectors": 1,
                                                  "gaia_weighted_variance_fraction_min": .4,
                                                  "gaia_min_accepted_epochs": 20,
                                                  "gaia_min_baseline_days": 100},
                      "cases": results}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
