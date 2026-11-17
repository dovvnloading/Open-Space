"""Independent, offline referee checks of the deposited Planet Nine campaign.

Reads the original joined trial rows, detector receipts, raw SCANPI members, and
catalogue response rather than trusting only the descriptive READMEs. Native
survey FITS files are *not* deposited; their original hashes remain receipts,
not newly authenticated downloads.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import tarfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
Z = 1.959963984540054


def load(root, relative):
    return json.loads((root / relative).read_text())


def records(root, relative):
    p = root / relative
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def same(actual, expected, atol=1e-8):
    assert math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=atol), (actual, expected)


def wilson(k, n):
    p = k / n
    denom = 1 + Z * Z / n
    center = (p + Z * Z / (2 * n)) / denom
    delta = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / denom
    return (center - delta, center + delta)


def deep_compare(current, baseline, label="", tol=5e-4):
    if isinstance(current, dict):
        assert isinstance(baseline, dict) and current.keys() == baseline.keys(), label
        for key in current:
            deep_compare(current[key], baseline[key], label + "/" + str(key), tol)
    elif isinstance(current, list):
        assert isinstance(baseline, list) and len(current) == len(baseline), label
        for j, (a, b) in enumerate(zip(current, baseline)):
            deep_compare(a, b, label + "/" + str(j), tol)
    elif isinstance(current, (int, float)) and not isinstance(current, bool):
        same(current, baseline, tol)
    else:
        assert current == baseline, label


def verify(root=ROOT):
    summary = load(root, "orbits/euclid/paired_summary.json")
    src = root / "archives/planet_draws.csv.gz"
    assert hashlib.sha256(src.read_bytes()).hexdigest() == summary["des_input_sha256"]
    des = {(r["experiment"], int(r["draw"])): r for r in records(root, "archives/planet_draws.csv.gz")}
    paired = records(root, "orbits/euclid/paired_trials.csv.gz")
    assert len(des) == len(paired) == 800
    assert {(r["experiment"], int(r["draw"])) for r in paired} == set(des)
    expected = {"published_centroids": (397, 379, 362, 344),
                "centroid_stress": (388, 371, 271, 254)}
    model = "mass5_R2.71_p0.10"
    for experiment, cases in summary["cases"].items():
        rs = [r for r in paired if r["experiment"] == experiment]
        assert len(rs) == 400
        for name, vals in cases.items():
            pairs = []
            for r in rs:
                a, b = int(r["DES_found_" + name]), int(r["Euclid_found_" + name])
                assert a == int(des[(experiment, int(r["draw"]))]["hit_" + name])
                e = int(r["either_found_" + name])
                assert e == (a | b)
                pairs.append((a, b, e))
            counts = dict(n=400, des=sum(x[0] for x in pairs),
                          euclid=sum(x[1] for x in pairs),
                          both=sum(x[0] & x[1] for x in pairs),
                          either=sum(x[2] for x in pairs),
                          neither=400 - sum(x[2] for x in pairs))
            for key, value in counts.items():
                assert vals[key] == value, (experiment, name, key)
            for actual, computed in zip(vals["wilson95"], wilson(counts["either"], 400)):
                same(actual, computed, 1e-12)
        v = cases[model]
        assert (v["either"], v["des"], v["euclid"], v["both"]) == expected[experiment]
    misses = {(r["experiment"], int(r["draw"])) for r in paired
              if not int(r["either_found_" + model])}
    blue = records(root, "n395_cal/paired_misses.csv")
    assert {(r["experiment"], int(r["draw"])) for r in blue} == misses
    assert Counter(r["experiment"] for r in blue) == {"published_centroids": 3, "centroid_stress": 12}
    usable = [r for r in blue if int(r["eligible"])]
    assert len(usable) == 7
    assert Counter(r["experiment"] for r in usable) == {"published_centroids": 1, "centroid_stress": 6}
    assert not any(int(r["recovered"] or 0) for r in blue)
    same(max(float(r["real_snr"]) for r in usable), 1.0172858238220215, 1e-9)
    assert sum(int(r["recovery_over_64"] or 0) for r in usable) == 3
    n395 = load(root, "n395_cal/injection_summary.json")
    assert sum(v["valid"] for v in n395["original_data_search"]) == 1466
    assert sum(v["above_5p5"] for v in n395["original_data_search"]) == 0
    bymag = {float(r["N395_AB"]): r for r in n395["by_magnitude"]}
    assert [(bymag[m]["eligible"], bymag[m]["recovered"]) for m in (23.25, 23.5)] == [(129, 126), (129, 91)]

    hist = load(root, "historical/scanpi/track_summary.json")
    july, late = hist["3_July_sign+1"], hist["3_Aug_Sep_sign+1"]
    same(july["fixed_Jy"], .2809686749250723)
    same(late["fixed_Jy"], .23105131787009972)
    same(late["galaxy_Jy"], .22138085322040768)
    same(late["moving_min_med_max_Jy"][0], -.08120534359233768)
    same(late["moving_min_med_max_Jy"][2], -.012774299282740341)
    assert july["gated_models"] == late["gated_models"] == 711
    all_scan = load(root, "historical/scanpi/all_scan_tracks.json")
    assert all_scan["grid_size"] == 3617
    assert all_scan["grid_meeting_flux_threshold_in_both_groups"]["0.05"] == 0
    same(all_scan["fixed_July_late_Jy"][0], july["fixed_Jy"])
    same(all_scan["fixed_July_late_Jy"][1], late["fixed_Jy"])
    # Confirm source-level SCANPI receipts by hashing the actual member bytes.
    receipts = load(root, "historical/scanpi/receipts.json")
    checked_profiles = 0
    with tarfile.open(root / "historical/scanpi/profiles.tar.gz") as archive:
        for r in receipts:
            if not r["name"].startswith(("b3/", "b4/")):
                continue
            if r.get("status") != 200:
                continue
            buf = archive.extractfile(r["name"]).read()
            assert len(buf) == r["length"] and hashlib.sha256(buf).hexdigest() == r["sha256"], r["name"]
            checked_profiles += 1
    assert checked_profiles >= 20, checked_profiles
    orient = records(root, "historical/scanpi/sign_check.csv")
    assert len(orient) == 5
    assert all(abs(float(r["predicted_galaxy_along_arcmin"]) - float(r["measured_galaxy_along_arcmin"])) < .011 for r in orient)
    wise = next(line.split() for line in (root / "historical/iras_allwise.txt").read_text().splitlines()
                if line.lstrip().startswith("J022301.05-483040.5"))
    ra0, dec0 = load(root, "orbits/orbit_input.json")["iras_radec_deg"]
    ra, dec = float(wise[1]), float(wise[2])
    ra0, dec0, ra, dec = [math.radians(x) for x in (ra0, dec0, ra, dec)]
    h = math.sin((dec-dec0)/2)**2 + math.cos(dec0)*math.cos(dec)*math.sin((ra-ra0)/2)**2
    separation = math.degrees(2*math.asin(math.sqrt(h)))*3600
    same(separation, 32.856107, .003)
    same(wise[8], 6.899, 1e-5)
    same(wise[9], 15.6, 1e-5)
    assert wise[10] == "0000"

    sp = load(root, "spherex/combined.json")
    stress = load(root, "spherex/centroid_stress/combined.json")
    assert sp["exposures"] == stress["exposures"] == 572
    assert len(sp["detectors"]) == len(stress["detectors"]) == 6
    assert (sp["grid_size"], sp["grid_valid"], stress["grid_size"], stress["grid_valid"]) == (3617, 2551, 15833, 10464)
    assert sp["grid_passes"] == stress["grid_passes"] == sp["control_passes"] == stress["control_passes"] == 0
    # One successfully cached D3 file has length/hash and UTC but lacks an HTTP status field.
    cached_missing_status = 0
    used_images = 0
    for index in range(1, 7):
        directory = "spherex/" + ("" if index == 1 else "D" + str(index) + "/")
        sr = load(root, directory + "summary.json")
        rs = load(root, directory + "images.json")
        usable = [r for r in rs if r.get("sha256") and r.get("length")]
        assert len(usable) == sr["exposures"], (index, len(usable), sr["exposures"])
        for r in usable:
            assert re.fullmatch("[0-9a-f]{64}", r["sha256"]) and r["length"] > 0
            if r.get("status") is None:
                cached_missing_status += 1
            else:
                assert r["status"] == 200
        used_images += len(usable)
    assert used_images == 572 and cached_missing_status == 1

    exposures = records(root, "archives/dedicated_2026.csv")
    assert len(exposures) == 12
    byrelease = Counter((r["dateobs_center"][:10], r["release_date"]) for r in exposures)
    assert byrelease == {("2026-09-13", "2028-03-05"): 4,
                         ("2026-10-28", "2028-04-19"): 4,
                         ("2026-11-15", "2028-05-07"): 4}
    assert all(re.fullmatch("[0-9a-f]{32}", r["md5sum"]) for r in exposures)
    # The following products were not deposited. Check receipts' syntax, not their bytes.
    raw_inputs = load(root, "historical/inputs.json")
    assert all(r.get("url", "").startswith("https://") and r["length"] > 0
               and re.fullmatch("[0-9a-f]{64}", r["sha256"]) for r in raw_inputs)
    for sub in ("orbits", "wise_vet", "act_vis"):
        for line in (root / sub / "SHA256SUMS").read_text().splitlines():
            digest, fname = line.split(maxsplit=1)
            target = root / sub / fname.strip().lstrip("*")
            assert target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == digest, target

    assert "2006-12-26T21:03:07.5" == load(root, "orbits/orbit_input.json")["akari_midpoint_utc"]
    assert "15:59:18.072" in (root / "orbits/TIMING.md").read_text()
    assert "17:38:34.553" in (root / "orbits/TIMING.md").read_text()
    return dict(matched_trials=800, dark_published_recovered=397,
                dark_stress_recovered=388, SCANPI_raw_profiles_hashed=checked_profiles,
                galaxy_separation_arcsec=separation, usable_N395_misses=7,
                spherex_exposures=used_images, spherex_cached_receipts_missing_http_status=1,
                spherex_stress_eligible=stress["grid_valid"], spherex_stress_draws=stress["grid_size"],
                embargoed_DECam_frames=len(exposures),
                upstream_native_images_not_in_repo=True)


def compare_replay(root):
    targets = ["historical/scanpi/track_summary.json",
               "historical/scanpi/all_scan_tracks.json",
               "spherex/combined.json", "spherex/centroid_stress/combined.json",
               "joint/matched_geometry.json"]
    for path in targets:
        deep_compare(load(root, path), load(ROOT, path), path)
    for path in ["historical/scanpi/sign_check.csv", "joint/matched_recovery.csv"]:
        a, b = records(root, path), records(ROOT, path)
        assert len(a) == len(b), path
        for ar, br in zip(a, b):
            assert ar.keys() == br.keys()
            for key in ar:
                try:
                    same(ar[key], br[key], .01 if key.startswith("measured_") else 1e-5)
                except ValueError:
                    assert ar[key] == br[key], (path, key)
    return targets + ["historical/scanpi/sign_check.csv", "joint/matched_recovery.csv"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--replayed", type=Path, default=None)
    args = parser.parse_args()
    results = verify()
    if args.replayed:
        results["replayed_outputs_compared"] = compare_replay(args.replayed)
    print(json.dumps(results, indent=2))
