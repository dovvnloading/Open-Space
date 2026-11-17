"""Independent, offline PR #142 checks on archived tables and query receipts.

These tests recompute statistics from *deposited records*. They do not rerun
pixel-level FITS reductions or independently authenticate remote FITS SHA-256s.
Use the individual science scripts plus their FITS manifests for that step.
"""
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import re
from statistics import median
import unittest

ROOT = Path(__file__).resolve().parents[1] / "results" / "jwst-nearby-companions"
SHA = re.compile(r"^[0-9a-f]{64}$")


def rows(path):
    with (ROOT / path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def jload(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def first_qualifying(records, level, recovered, injected, threshold):
    """First grid point meeting the fractional recovery cut."""
    counts = defaultdict(lambda: [0, 0])
    for record in records:
        key = float(record[level])
        counts[key][0] += int(record[recovered])
        counts[key][1] += int(record[injected]) if injected else 1
    eligible = [
        k for k, (hits, total) in counts.items()
        if total and hits / total >= threshold
    ]
    return min(eligible) if eligible else None


class Referee(unittest.TestCase):

    def test_frozen_inventory(self):
        v = rows("inventory/archive/visits.csv")
        n = rows("inventory/archive/nearby-fields.csv")
        s = rows("inventory/archive/target-search-status.csv")
        self.assertEqual(len(v), 2653)
        self.assertEqual(sum(int(x["filesets"]) for x in v), 9056)
        self.assertEqual(len({x["program"] for x in v}), 240)
        self.assertEqual(len(n), 792)
        self.assertEqual(len({x["program"] for x in n}), 70)
        self.assertEqual(len({x["gcns_edr3_source_id"] for x in n}), 331)
        self.assertEqual(len(s), 331)
        self.assertEqual(sum(int(x["title_selected_papers"]) > 0 for x in s), 214)
        self.assertEqual(sum(int(x["title_selected_JWST_papers"]) > 0 for x in s), 83)

    def test_l788_pan_starrs_from_individual_epochs(self):
        p = rows("vetting/ps1_motion.csv")
        self.assertEqual(len(p), 14)
        self.assertEqual(len({x["date"] for x in p}), 9)
        self.assertAlmostEqual(median(float(x["fixed_resid_arcsec"]) for x in p),
                               .252157021625829, places=10)
        self.assertAlmostEqual(median(float(x["comoving_resid_arcsec"]) for x in p),
                               11.133507395064719, places=9)

    def test_g29_spitzer_individual_measurements(self):
        x = [p for p in jload("g29/spitzer/measurements.json") if p["band"] == "I2"]
        self.assertEqual(len(x), 2)
        x.sort(key=lambda p: p["date"])
        self.assertEqual([p["gaia_ties"] for p in x], [24, 24])
        for p, expected, motion in zip(x, [0.06970585598292939, 0.05526846294444105],
                                       [10.003018666251965, 7.834426660756949]):
            self.assertAlmostEqual(p["fixed_centroid_distance_arcsec"], expected, places=10)
            self.assertAlmostEqual(p["source_offset_fixed_vs_comoving_arcsec"], motion, places=9)

    def test_gj205_gaia_and_source_dispositions(self):
        f = jload("vetting/field_identifications.json")["hip"]
        measured, predicted = f["relative_arcsec"], f["gaia_relative_pred_arcsec"]
        self.assertAlmostEqual(math.hypot(*measured), 5.43647698422425, places=8)
        self.assertAlmostEqual(math.dist(measured, predicted), .05102884059501458, places=10)
        self.assertAlmostEqual(float(f["gaia_host"]["parallax"]), 175.31313655111924)
        self.assertAlmostEqual(float(f["gaia_field_star"]["parallax"]), 1.140226810130067)
        ross = [r for r in rows("nircam/nircam_trigger_checks.csv")
                if r["field"] == "ross619"]
        self.assertEqual(len(ross), 2)
        self.assertTrue(all(float(r["snr"]) > 5.7 for r in ross))
        self.assertTrue(all(float(r["archival_rdi_snr"]) < 3.1
                            and r["archival_above_cut"] == "0" for r in ross))
        wd = rows("miri/measured/sources.csv")
        self.assertEqual(len(wd), 16)
        statuses = Counter((r["target"], r["status"]) for r in wd)
        self.assertEqual(statuses["WD1202-232", "stationary match"], 4)
        self.assertEqual(statuses["WD1202-232", "unmatched"], 7)
        self.assertEqual(statuses["WD2105-82", "stationary match"], 5)
        self.assertFalse(any(r["status"] == "comoving" for r in wd))

    def test_nircam_injection_recovery(self):
        totals = 0
        for name in ("g99_49", "l788_34", "ross619", "v374peg", "mcc124"):
            c = rows("nircam/" + name + "_contrast.csv")
            inj = rows("nircam/" + name + "_injections.csv")
            self.assertEqual(len(c), 6)
            self.assertEqual(len(inj), 2376)
            totals += len(inj)
            for cur in c:
                same = [r for r in inj if r["radius_arcsec"] == cur["radius_arcsec"]]
                self.assertEqual(len(same), 36 * 11)
                by_level = defaultdict(list)
                for r in same:
                    by_level[float(r["flux_uJy"])].append(int(r["recovered"]))
                lv = sorted(by_level)
                self.assertEqual(len(lv), 11)
                candidates = [f for i, f in enumerate(lv)
                              if all(sum(by_level[g]) >= 35 for g in lv[i:])]
                self.assertTrue(candidates, (name, cur["radius_arcsec"]))
                self.assertAlmostEqual(min(candidates), float(cur["f95_uJy"]), places=10)
        self.assertEqual(totals, 11880)
        g = rows("nircam/g99_49_contrast.csv")
        self.assertAlmostEqual(min(float(x["contrast95"]) for x in g), 1.7732407393360196e-6)
        self.assertAlmostEqual(max(float(x["contrast95"]) for x in g), 2.695220583511479e-6)

    def test_gj205_miri_mask_recovery(self):
        for band in ("f1065c", "f1550c"):
            cs = rows("gj205/" + band + "_miri_contrast.csv")
            ins = rows("gj205/" + band + "_miri_injections.csv")
            self.assertEqual(len(cs), 6)
            self.assertEqual(len(ins), 6 * 48 * 11)
            for c in cs:
                selection = [r for r in ins
                             if r["radius_arcsec"] == c["radius_arcsec"]
                             and r["eligible"] == "1"]
                self.assertEqual(len(selection), int(c["eligible_angles"]) * 11)
                by_level = defaultdict(list)
                for r in selection:
                    by_level[float(r["flux_uJy"])].append(int(r["recovered"]))
                lv = sorted(by_level)
                self.assertEqual(len(lv), 11)
                candidates = [f for i, f in enumerate(lv)
                              if all(sum(by_level[g]) / len(by_level[g]) >= .95
                                     for g in lv[i:])]
                self.assertTrue(candidates)
                self.assertAlmostEqual(min(candidates), float(c["f95_uJy"]), places=8)

    def test_other_miri_recovery(self):
        for folder, curve, trials, key, level, count, nrows in (
            ("g29/out", "contrast_recovery.csv", "injection_trials.csv",
             ("filter", "radius_arcsec"), "contrast", 1602, 48),
            ("ross", "contrast_recovery.csv", "injection_trials.csv",
             ("obs", "radius_arcsec"), "injected_contrast", 1220, 60)):
            cs = rows(folder + "/" + curve)
            t = rows(folder + "/" + trials)
            self.assertEqual(len(t), count)
            self.assertEqual(len(cs), nrows)
            groups = defaultdict(list)
            for r in cs:
                groups[tuple(r[x] for x in key)].append(r)
                self.assertGreater(int(r["injected"] if folder == "ross"
                                       else r["tried"]), 0)
            for group, matches in groups.items():
                reco = "recovered"
                den = "injected" if folder == "ross" else "tried"
                got = first_qualifying(matches, level, reco, den, .90)
                self.assertIsNotNone(got, group)
                self.assertGreater(got, 0)
        wd = rows("miri/measured/contrast.csv")
        self.assertEqual(len(rows("miri/measured/injections.csv")), 1920)
        self.assertEqual(len(wd), 48)
        expected = {
            ("WD1202-232", 5): .01, ("WD1202-232", 8): .01,
            ("WD1202-232", 11): .04, ("WD2105-82", 5): .04,
            ("WD2105-82", 8): .02, ("WD2105-82", 11): .02,
        }
        for (star, rad), truth in expected.items():
            group = [r for r in wd if r["target"] == star
                     and float(r["radius_arcsec"]) == rad]
            self.assertAlmostEqual(first_qualifying(group, "companion_host_ratio",
                                                     "n_recovered", "n_injected", .90),
                                   truth, places=8)

    def test_mast_query_receipt_structure(self):
        for key, expected in (("coron", 5299), ("stellar_imaging", 4192),
                              ("nearby_programs", 1008)):
            d = jload("inventory/archive/" + key + "-receipt.json")
            self.assertEqual(d["returned"], expected)
            self.assertEqual(d["public_at_cutoff"], expected)
            self.assertIn("mast.stsci.edu", d["endpoint"])
            self.assertEqual(sum(int(p["returned"]) for p in d["pages"]), expected)
            self.assertTrue(all(SHA.fullmatch(p["sha256_response"]) for p in d["pages"]))
            self.assertTrue(all("publicReleaseDate" in str(p["request"]) for p in d["pages"]))
        # URI, PUBLIC rights and SHA-256 syntax only; binary contents must be
        # independently fetched and hashed before treating these as validated.
        inputs = jload("g29/inputs.json")
        self.assertEqual(sum(r["proposal_id"] == "9077" for r in inputs), 2)
        self.assertTrue(all(SHA.fullmatch(r["sha256"]) and
                            r["dataRights"] == "PUBLIC" for r in inputs))
        for p in rows("nircam/go8826_inputs.csv"):
            self.assertTrue(SHA.fullmatch(p["sha256"]))
            self.assertEqual(p["rights"], "PUBLIC")


if __name__ == "__main__":
    unittest.main()
