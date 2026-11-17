"""Independent, offline referee checks for the bounded Euclid Q1 dark-lens null.

Compute the printed claims directly from deposited tabular inputs; do not import
any analysis code from results/euclid_q1_dark_lenses.
"""
import collections
import csv
import hashlib
import json
import math
import statistics
import unittest
import zipfile
from pathlib import Path

R = Path(__file__).resolve().parents[1] / "results" / "euclid_q1_dark_lenses"
B = "102158888_2659753895653871052"


def table(name):
    with (R / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def close(actual, expected, tol=0.005):
    if not math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=tol):
        raise AssertionError(f"Expected {expected}, got {actual}")


class EuclidReferee(unittest.TestCase):
    def test_counts_and_minimum_source_subtracted_flux(self):
        population = table("catalog/slde_screen.csv")
        self.assertEqual(len(population), 500)
        self.assertEqual(len({r["id_str"] for r in population}), 500)
        self.assertEqual(collections.Counter(r["model_state"] for r in population),
                         {"positive_radius": 312, "no_row": 174, "missing_radius": 12,
                          "nonpositive_radius": 2})
        self.assertEqual(collections.Counter(r["grade"] for r in population),
                         {"A": 253, "B": 247})
        modeled = {r["id_str"] for r in population if r["model_state"] == "positive_radius"}
        raw = table("imaging/apertures.csv")
        cores = [r for r in raw if r["id_str"] in modeled and r["band"] == "VIS"
                 and r["aperture"] == "core" and r["coverage_ok"] == "True"]
        self.assertEqual(len(cores), 312)
        self.assertEqual({r["id_str"] for r in cores}, modeled)
        lowest = min(float(r["source_subtracted_snr"]) for r in cores)
        close(lowest, 19.74, 0.005)
        self.assertTrue(all(float(r["source_subtracted_snr"]) > 0 for r in cores))
        for band, expected in [("VIS", 26.12), ("NIR_Y", 24.10),
                               ("NIR_J", 24.04), ("NIR_H", 23.88)]:
            depths = [float(r["noise_only_5sigma_mag"]) for r in raw
                      if r["id_str"] in modeled and r["band"] == band
                      and r["aperture"] == "core" and r["coverage_ok"] == "True"
                      and r["noise_only_5sigma_mag"]]
            close(statistics.median(depths), expected, 0.01)
        self.assertEqual(len({r["id_str"] for r in population} - modeled), 188)

    def test_blind_area_and_trigger_receipts(self):
        strips = table("arc_search/strips.csv")
        self.assertEqual(len(strips), 32)
        self.assertEqual(len({(r["tile"], r["row"]) for r in strips}), 32)
        self.assertEqual(len({r["tile"] for r in strips}), 4)
        eligible = sum(int(r["eligible_binned_centers"]) for r in strips)
        self.assertEqual(eligible, 27013228)
        for row in strips:
            self.assertEqual(int(row["last_byte"]) - int(row["first_byte"]) + 1,
                             int(row["height"]) * 19200 * 4)
        close(eligible * (0.2 / 3600) ** 2, 0.083374, 0.0000005)
        triggers = table("arc_search/vis_triggers.csv")
        self.assertEqual(len(triggers), 1258)
        self.assertEqual(sum(int(r["pass_vis"]) for r in triggers), 6)
        with zipfile.ZipFile(R / "arc_search/frozen_stages.zip") as z:
            self.assertIsNone(z.testzip())
            self.assertEqual(len(z.namelist()), 64)
            for row in strips:
                key = row["tile"] + "_" + row["row"] + ".json"
                run = json.loads(z.read("runs/" + key))
                profile = json.loads(z.read("profiles/" + key))
                self.assertEqual(run["sha256"], row["sha256"])
                self.assertEqual(len(run["triggers"]), int(row["n_triggers"]))
                self.assertEqual(len(profile["profiles"]), int(row["n_triggers"]))
                self.assertEqual(sum(int(r["pass_vis"]) for r in profile["profiles"]),
                                 int(row["n_profiles_passing"]))
        searched = json.loads((R / "arc_search/unmodelled.json").read_text())
        fallback = json.loads((R / "arc_search/fallback.json").read_text())
        self.assertEqual(len(fallback), 13)
        self.assertTrue(all(r["status"] == "examined" for r in fallback))
        examined = {r["id_str"] for r in searched + fallback if r["status"] == "examined"}
        modeled = {r["id_str"] for r in table("catalog/slde_screen.csv")
                   if r["model_state"] == "positive_radius"}
        self.assertEqual(examined, {r["id_str"] for r in table("catalog/slde_screen.csv")} - modeled)
        self.assertEqual(sum(int(r.get("pass_vis", 0)) for r in searched + fallback), 1)

    def test_injection_response_not_an_abundance_limit(self):
        rows = table("arc_search/lens_injections.csv")
        self.assertEqual(len(rows), 1491)
        self.assertTrue(all(r["preexisting"] == "0" for r in rows))
        trial = collections.Counter((r["source_offset_over_theta"],
                                      r["intrinsic_peak_sigma"]) for r in rows)
        recovered = collections.Counter((r["source_offset_over_theta"],
                                          r["intrinsic_peak_sigma"])
                                         for r in rows if r["recovered"] == "1")
        expected = {("0.02", "3"): (127, 248),
                    ("0.02", "8"): (213, 249),
                    ("0.18", "3"): (28, 248),
                    ("0.18", "8"): (150, 249),
                    ("0.45", "3"): (0, 248),
                    ("0.45", "8"): (0, 249)}
        for k, (yes, total) in expected.items():
            self.assertEqual((recovered[k], trial[k]), (yes, total))
        self.assertEqual(sum(recovered.values()), 518)
        self.assertEqual(len(table("arc_search/validation/conditional_limits.csv")), 700)

    def test_secondary_redshift_mass_to_light_and_counterimage(self):
        population = {r["id_str"]: r for r in table("catalog/slde_screen.csv")}
        lead = population[B]
        close(lead["right_ascension"], 265.975390, .000001)
        close(lead["declination"], 65.387105, .000001)
        close(lead["theta_eff_arcsec_median_pdf"], .93605511, .000001)
        rows = {(r["id_str"], r["band"], r["z_l_model"], r["zs"]): r
                for r in table("inference/conditional_ML_posterior.csv")}
        for model, zsource, med, low, high in [
            ("Euclid_PHZ", "2", 105.2, 89.0, 123.1),
            ("Legacy_DR9", "2", 54.4, 35.1, 76.3),
            ("PHZ_PP_median_fixed", "3", 10.7, 10.1, 11.5),
            ("PHZ_PP_median_fixed", "4", 5.02, 4.73, 5.37)]:
            row = rows[B, "h", model, zsource]
            tol = 0.055 if med >= 10 else .0051
            close(row["ML_p50"], med, tol)
            close(row["ML_p16"], low, .055 if low >= 10 else .0051)
            close(row["ML_p84"], high, .055 if high >= 10 else .0051)
            self.assertEqual(int(row["n"]), 48000)
        counter = next(r for r in table("imaging/secondary_counterimage.csv")
                       if r["band"] == "VIS" and r["light_model"] == "MGE"
                       and r["aperture_arcsec"] == "0.3")
        close(counter["predicted_counter_snr"], 11.46, .005)
        close(counter["counter_snr"], -1.19, .005)
        close(float(counter["flux_ratio_sis"]) *
              float(counter["primary_counts"]) /
              (float(counter["counter_counts"]) / float(counter["counter_snr"])),
              counter["predicted_counter_snr"], 0.00005)
        self.assertEqual(next(r for r in table("imaging/secondary_models.csv")
                              if r["model"] == "SIS+shear")["usable_pixels"], "1396")
        physical = json.loads((R / "inference/lead_physical_query.json").read_text())
        self.assertEqual(hashlib.sha256(physical["raw_csv"].encode("utf-8")).hexdigest(),
                         physical["response_sha256"])
        self.assertIn("2.4572818279266357", physical["raw_csv"])
        self.assertIn("0.2081272304058075", physical["raw_csv"])
        for file in ("inference/lead_spectra_query.json", "inference/lead_physical_query.json"):
            source = json.loads((R / file).read_text())
            self.assertEqual(hashlib.sha256(source["raw_csv"].encode()).hexdigest(),
                             source["response_sha256"], file)
        for rec in json.loads((R / "inference/cone_responses.json").read_text()):
            self.assertEqual(hashlib.sha256(rec["raw_csv"].encode()).hexdigest(),
                             rec["sha256"], rec.get("id", rec.get("survey", "cone")))

    def test_competing_models_and_source_receipt_schema(self):
        rows = table("arc_search/model/cross_validation.csv")
        nominal = [r for r in rows if r["center_bound_arcsec"] == "0.35"]
        self.assertEqual(len(nominal), 7)
        self.assertEqual(sum(float(r["mean_square_knots"]) <
                             float(r["mean_square_lens"]) for r in nominal), 5)
        self.assertEqual(len(table("arc_search/model/cross_validation_folds.csv")), 45)
        native = table("arc_search/validation/native_inputs.csv")
        self.assertEqual(len(native), 84)
        self.assertEqual(len({r["sha256"] for r in native}), 84)
        self.assertEqual(len(table("imaging/input_members.csv")), 1288)
        peers = next(r for r in table("inference/redshift_matched_peers.csv")
                     if r["id_str"] == B and r["z_l_model"] == "Legacy_DR9")
        self.assertEqual(peers["n_peers"], "25")
        close(peers["peer_ML_p50"], 2.64, .005)
        close(peers["peer_ML_max"], 18.32, .005)
        close(peers["A_H_to_peer_median_mag"], 3.29, .005)
        lens = next(r for r in table("imaging/secondary_models.csv")
                    if r["model"] == "SIS+shear")
        knot = table("imaging/secondary_blend.csv")[0]
        close(lens["radius_arcsec"], .480, .0005)
        close(lens["light_mass_offset_arcsec"], .766, .0005)
        close(lens["residual_diagonal_chi2"], 915.9, .05)
        close(knot["diagonal_chi2"], 920.1, .05)
        self.assertEqual((lens["usable_pixels"], knot["usable_pixels"]), ("1396", "1396"))
        for r in native + table("imaging/input_members.csv"):
            self.assertEqual(len(r["sha256"]), 64)
            int(r["sha256"], 16)
            self.assertGreater(int(r["bytes"]), 0)
        # The deposited SHA256SUMS checks committed bytes. External input receipts
        # are metadata: checking their syntax does not independently hash remote bytes.

if __name__ == "__main__":
    unittest.main(verbosity=2)
