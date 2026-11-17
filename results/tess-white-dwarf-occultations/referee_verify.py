"""Independent, offline referee check of frozen TESS WD occultation claims.

Reads committed inputs, not report prose. The specialised audit scripts separately
replay random draws, period endpoint matching, and the finite-trial model.
"""
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import beta

R = Path(__file__).resolve().parent


def read(name):
    path = R / name
    with (gzip.open(path, "rt") if name.endswith(".gz") else open(path)) as stream:
        return json.load(stream)


def table(name):
    path = R / name
    with (gzip.open(path, "rt") if name.endswith(".gz") else open(path)) as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256((R / path).read_bytes()).hexdigest()


def check():
    selection = read("extension/extension_selection.json")
    parents = table("extension/extension_parent.csv")
    assert len(parents) == 157
    assert sha("extension/extension_parent.csv") == selection["roster_sha256"]
    assert Counter(x["selection"] for x in parents) == selection["decisions"]
    assert selection["decisions"]["request latest visit"] == 140
    measurements = read("extension/extension_measurements.json")
    assert len(measurements) == 140
    status = Counter(m["status"] for m in measurements)
    assert status == {
        "included": 39,
        "excluded: Gaia crowding cut": 89,
        "excluded: nonpositive aperture baseline": 4,
        "no available 200-second visit": 8,
    }
    accepted = [m for m in measurements if m["status"] == "included"]
    assert len({m["source_id"] for m in accepted}) == 39
    assert all(m["min_target_fraction"] >= .8 and m["gaia_without_G"] == 0
               and 56 <= m["sector"] <= 108 and abs(m["cadence_s"]-200) < 1
               for m in accepted)
    assert sum(m["quality0_points"] for m in accepted) == 417350
    assert sum(m["search"]["usable_cadences"] for m in accepted) == 417173
    assert sum(len(m["search"]["periodic"]) for m in accepted) == 19
    assert sum(len(m["search"]["isolated"]) for m in accepted) == 85
    assert sum(bool(m["search"]["periodic"] or m["search"]["isolated"]) for m in accepted) == 20
    receipts = {r["file"]: r for r in read("extension/extension_receipts.json")}
    for m in accepted:
        entry = receipts[f'{m["source_id"]}-s{m["sector"]}.fits.json']
        assert entry["sha256"] == m["fits_sha256"] and entry["http_status"] == 200
    assert all(len(r.get("sha256", "")) == 64 for r in receipts.values()
               if r.get("status") == "ok")

    cells = read("extension/extension_limits.json")
    stars = table("extension/extension_recovery_stars.csv")
    assert len(cells) == 4 and len(stars) == 39 * 4
    for c in cells:
        rows = [s for s in stars
                if float(s["period_days"]) == c["period_days"]
                and float(s["radius_earth"]) == c["radius_earth"]]
        assert len(rows) == 39 and c["targets"] == 39
        assert sum(int(s["possible_host"]) for s in rows) == c["possible_hosts"] == 0
        assert sum(int(s["recovered"]) for s in rows) == c["recovered"]
        assert c["trials"] == 1950
        bounds = [beta.ppf(.025 / 39, int(s["recovered"]), 51-int(s["recovered"]))
                  if int(s["recovered"]) else 0 for s in rows]
        func = lambda f: math.prod(1-f*q for q in bounds) - .025
        assert func(1) < 0
        assert math.isclose(brentq(func, 0, 1), c["transiting_host_upper95"],
                            rel_tol=1e-10)
        assert c["intrinsic_host_upper95"] is None
    target_cell = next(c for c in cells if c["radius_earth"] == 11.2
                       and c["period_days"] == 1.40794)
    assert target_cell["recovered"] == 1896
    assert abs(100*target_cell["transiting_host_upper95"]-11.130736014663589) < 1e-8

    population = read("extension/population_limit.json")
    pstars = read("extension/population_stars.json")
    assert len(pstars) == population["targets"] == 42
    assert (population["period_min_days"], population["period_max_days"],
            population["radius_earth"]) == (.25, 3, 11.2)
    assert population["trials"] == 2100 and population["recovered"] == 2035
    possible = {s["target"] for s in pstars if s["possible_host"]}
    assert possible == set(population["possible_host_ids"]) == {
        "4672306015773211008", "5063539946887524864"}
    assert population["possible_hosts"] == len(possible) == 2
    periods = {m["source_id"]: [r["period_days"] for r in m["search"]["periodic"]]
               for m in accepted}
    for star in pstars:
        if star["target"] in possible:
            for p in star["compatible_periods_days"]:
                assert p in periods[star["target"]] and .25 <= p <= 3
    lows = [beta.ppf(.025/42, s["recovered"], 51-s["recovered"])
            if s["recovered"] else 0 for s in pstars]

    def tail(f):
        # Direct Poisson-binomial P(K<=2), not the survey's odds-polynomial code.
        a, b, c = 1., 0., 0.
        for q in lows:
            v = f*q
            a, b, c = a*(1-v), b*(1-v)+a*v, c*(1-v)+b*v
        return a+b+c

    assert math.isclose(brentq(lambda f: tail(f)-.025, 0, 1),
                        population["transiting_host_upper95"], rel_tol=1e-10)
    assert abs(100*population["transiting_host_upper95"] - 20.125634945490806) < 1e-8
    assert population["intrinsic_host_upper95"] is None
    assert math.isclose(population["geometry_weighted_stars"], 1.6780708988402504)

    control = read("controls/planet_target_pixel_results.json")
    assert len(control["products"]) == 2
    assert {p["cadence_s"] for p in control["products"]} == {20, 120}
    for p in control["products"]:
        assert p["tic"] == 267574918 and p["sector"] == 86
        assert len(p["sha256"]) == 64 and p["bytes"] > 1000000
        selected = p["apertures"]["flux_error_optimized"]
        hit = selected["matching_peak"]
        assert hit["rank"] == 1 and abs(hit["period_days"] - 1.407939217) < .0001
        assert hit["events_covered"] >= 10 and hit["score"] > 10
        assert selected["strongest_peak"]["period_days"] == hit["period_days"]
        assert p["apertures"]["radius_1p5"]["matching_peak"] is None

    recent = read("sector110/tic55997020/search.json")
    assert recent["quality0_finite"] == 2123 and recent["attempted"] == 2549
    assert len(recent["search"]["periodic"]) == 0
    assert len(recent["search"]["isolated"]) == 1
    hit = recent["search"]["isolated"][0]
    assert hit["duration_min"] == 4 and abs(hit["depth"]-.15646598677999404)<1e-12
    vet = read("sector110/tic55997020/image_vet.json")[0]
    assert vet["event"] == hit
    fixed = {a["radius"]: a for a in vet["target"]["fixed_apertures"]}
    for r, a in fixed.items():
        assert math.isclose(1-a["in_flux"]/a["shoulder_flux"], a["depth"], abs_tol=1e-12)
    assert abs(fixed[1.5]["depth"] - .20659731512474397) < 1e-12
    assert abs(fixed[3.5]["depth"] - .001455419175957854) < 1e-12
    assert vet["comparison"]["fixed_apertures"][0]["depth"] > .25
    assert abs(vet["target"]["centroid_change_pix"][0]) > .10
    assert abs(vet["comparison"]["centroid_change_pix"][0]) > .10

    screens = read("short_cadence/roster_screen.json.gz")
    used = [r for r in screens if r.get("status") == "searched"]
    assert len(screens) == 1693 and len(used) == 1106
    assert len({r["tic"] for r in used}) == 86
    assert sum(len(r["periodic"]) for r in used) == 0
    assert sum(len(r["isolated"]) for r in used) == 142

    optics = read("vetting/catalogue_checks.json.gz")
    assert len(optics) == 37
    assert sum(len(r["events"]) for r in optics) == 272
    optical_curves = read("vetting/optical_epochs.json.gz")
    assert sum(len(c["btjd"]) for c in optical_curves if c["band"] == "ztf") == 32719
    ground = read("ground_windows/coverage.json")
    assert ground["matching_boxes"] == ground["matching_images"] == 0
    assert ground["targets"] == 37 and ground["groups"] == 59
    assert sum(ground["window_counts"].values()) == 252
    assert ground["distinct_atlas_exposures"] == 336
    assert sha("ground_windows/responses.json.gz") == ground["responses_sha256"]
    assert sha("ground_windows/windows.json.gz") == ground["windows_sha256"]

    metrics = dict(
        parent=157, requested=140, crowding_excluded=89, included=39,
        quality0=417350, detector_cadences=417173, extension_periodic=19,
        extension_isolated=85, planet_injections="1896/1950",
        fixed_period_days=1.40794, fixed_radius_earth=11.2,
        fixed_period_transiting_upper95_pct=100*target_cell["transiting_host_upper95"],
        continuous_period_range_days=[.25, 3], continuous_possible_hosts=sorted(possible),
        continuous_transiting_upper95_pct=100*population["transiting_host_upper95"],
        control_pixel_cadences=[20,120], sector110_search_depth_pct=100*hit["depth"],
        sector110_fixed_1p5_depth_pct=100*fixed[1.5]["depth"],
        sector110_fixed_3p5_depth_pct=100*fixed[3.5]["depth"],
        fast_products=len(used), fast_stars=86, optical_hosts=len(optics),
        optical_trigger_records=272, ztf_matched_points=32719)
    print("REFEREE_VERIFY " + json.dumps(metrics, sort_keys=True))


if __name__ == "__main__":
    check()
