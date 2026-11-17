"""Independent offline referee replay of the deposited PR #129 scientific tables.

No archive or network access is needed. Scripts are run on a scratch copy so
the repository's frozen data and receipts are not modified by validation.
Raw FITS mosaic photometry remains a separate, network-dependent replay.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy.integrate import quad

SRC = Path(__file__).resolve().parents[1] / "results/jwst-repeat-imaging-transients"


def check_csv_equal(testcase, before, after):
    testcase.assertEqual(list(before.columns), list(after.columns))
    pd.testing.assert_frame_equal(
        before, after, check_dtype=False, check_exact=False,
        rtol=3e-5, atol=3e-5,
    )


class JWSTPR129Referee(unittest.TestCase):
    def test_reexecute_deposited_scripts_in_scratch_copy(self):
        with tempfile.TemporaryDirectory(prefix="jwst129-referee-") as temp:
            target = Path(temp) / "survey"
            shutil.copytree(SRC, target)
            jobs = [
                ("repeat-catalog-followup/cosmos-022-history", "summarize.py",
                 ["nightly.csv", "medium_changes.csv", "inventory_checksums.csv"]),
                ("repeat-catalog-followup/cosmos-022-validation", "emission_lines.py", ["muse_lines.csv"]),
                ("repeat-catalog-followup/cosmos-022-validation", "broad_hbeta.py", ["muse_broad_hbeta.csv"]),
                ("repeat-catalog-followup", "dispositions.py", ["dispositions.csv"]),
                ("repeat-catalog-followup", "candidate_photometry.py", ["cosmos022_host_magnitude.json"]),
                ("goods-s-oct2026", "blind_summary.py",
                 ["blind_recovery.csv", "blind_area.csv", "blind_pairs.csv"]),
                ("goods-s-oct2026", "nuclear.py", ["nuclear_candidates.csv"]),
                ("high-z-hosts", "join_bands.py", ["joint_contrasts.csv.gz",
                 "joint_injections.csv.gz", "joint_recovery.csv"]),
            ]
            for folder, script, outputs in jobs:
                cwd = target / folder
                before = {f: json.loads((cwd / f).read_text()) if f.endswith(".json")
                          else pd.read_csv(cwd / f) for f in outputs}
                cmd = [sys.executable, script]
                if script == "join_bands.py":
                    cmd.append(".")
                completed = subprocess.run(
                    cmd, cwd=cwd, capture_output=True, text=True, timeout=120,
                    env={**os.environ, "MPLBACKEND": "Agg"},
                )
                self.assertEqual(completed.returncode, 0,
                                 f"{folder}/{script}: {completed.stderr[-1800:]}")
                for f in outputs:
                    with self.subTest(script=script, output=f):
                        if f.endswith(".json"):
                            out = json.loads((cwd / f).read_text())
                            self.assertEqual(before[f].keys(), out.keys())
                            for k, v in before[f].items():
                                if isinstance(v, (int, float)):
                                    self.assertAlmostEqual(v, out[k], places=5, msg=k)
                                else:
                                    self.assertEqual(v, out[k])
                        else:
                            check_csv_equal(self, before[f], pd.read_csv(cwd / f))

    def test_catalog_census_and_cosmos_flux(self):
        root = SRC / "repeat-catalog-followup"
        inputs = pd.read_csv(root / "catalog_inputs.csv")
        counts = pd.read_csv(root / "screen_counts.csv")
        groups = pd.read_csv(root / "coincidences.csv")
        dispositions = pd.read_csv(root / "dispositions.csv")
        self.assertEqual((len(inputs), inputs["sources"].sum()), (69, 317286))
        self.assertEqual((len(counts), int(counts.eligible.sum()),
                          int(counts.triggers.sum())), (37, 47584, 463))
        self.assertEqual(len(groups), 42)
        self.assertEqual(set(groups.candidate), set(dispositions.candidate))
        self.assertEqual(dispositions.bands_passing_both_apertures.ge(2).sum(), 14)
        self.assertEqual(len(pd.read_csv(root / "pixel_apertures.csv")
                             [["candidate", "band"]].drop_duplicates()), 100)
        history = root / "cosmos-022-history"
        night = pd.read_csv(history / "nightly.csv")
        self.assertEqual((len(night), night.utc_night.nunique(),
                          night["filter"].nunique()), (105, 20, 20))
        for band in ["F140M", "F162M", "F182M", "F210M",
                     "F250M", "F300M", "F360M", "F460M"]:
            m = pd.read_csv(history / "medium_changes.csv", index_col=0)
            self.assertGreater(m.loc[band, "delta_njy"], 0)
        params = json.loads((root / "cosmos022_host_magnitude.json").read_text())
        self.assertAlmostEqual(params["delta_aperture_njy"], 236.179315, places=5)
        self.assertAlmostEqual(params["empirical_error_njy"], 18.2375653, places=5)
        self.assertAlmostEqual(params["spectroscopic_host_absolute_difference_ab"], -16.44927, places=4)
        self.assertEqual(params["host_id"], 953994)
        self.assertAlmostEqual(params["spec_z"], 0.5067, places=5)

    def test_native_exposures_and_spectrum(self):
        v = SRC / "repeat-catalog-followup/cosmos-022-validation"
        c = pd.read_csv(v / "native_cal.csv")
        self.assertEqual(len(c), 34)
        for band, old, new in [("f162m", 611, 843), ("f182m", 624, 824)]:
            a = c[(c.band == band) & c.mosaic.str.contains("-o004")]
            b = c[(c.band == band) & c.mosaic.str.contains("-o003")]
            self.assertEqual((len(a), len(b)), (9, 8))
            self.assertAlmostEqual(a.flux_r16_njy.median(), old, delta=1)
            self.assertAlmostEqual(b.flux_r16_njy.median(), new, delta=1)
            self.assertGreater(b.flux_r16_njy.min(), a.flux_r16_njy.max())
            a = a[a.DO_NOT_USE == 0]
            b = b[b.DO_NOT_USE == 0]
            self.assertEqual((len(a), len(b)), (6, 4))
            self.assertGreater(b.flux_r16_njy.min(), a.flux_r16_njy.max())
        for bit in ["SATURATED", "DROPOUT", "PERSISTENCE", "OUTLIER"]:
            self.assertEqual(int(c[bit].sum()), 0, bit)
        lines = pd.read_csv(v / "muse_lines.csv").set_index("line")
        self.assertAlmostEqual(lines.loc["OIII5007", "flux_1e20"] /
                               lines.loc["Hbeta", "flux_1e20"], 2.91, delta=.06)
        self.assertAlmostEqual(lines.loc["OII3727", "flux_1e20"] /
                               lines.loc["Hbeta", "flux_1e20"], 1.73, delta=.06)

    def test_goods_candidate_and_blind_depth(self):
        g = SRC / "goods-s-oct2026"
        stats = pd.read_csv(g / "catalog_stats.csv")
        self.assertEqual((len(stats), int(stats.mutual.sum())), (12, 20568))
        params = pd.read_csv(g / "nuclear_candidates.csv").set_index("site").loc[185]
        self.assertEqual(int(params.dr2_id), 103607)
        self.assertAlmostEqual(params.zphot, 2.47)
        self.assertAlmostEqual(params.zphot_68lo, 1.907)
        self.assertAlmostEqual(params.zphot_68hi, 2.514)
        self.assertAlmostEqual(params.difference_f277w_ap_njy, 22.2483622, places=5)
        self.assertAlmostEqual(params.diff_local_scatter_njy, 1.15446758, places=5)
        self.assertAlmostEqual(params.abs_delta_ab_ap, -17.18484, places=4)
        self.assertEqual((int(params.old_dither_count), int(params.new_dither_count)), (2, 2))
        self.assertAlmostEqual(params.old_dither_median_njy, 50.6104, places=3)
        self.assertAlmostEqual(params.new_dither_median_njy, 71.4944, places=3)
        area = pd.read_csv(g / "blind_area.csv").set_index("filter")
        self.assertAlmostEqual(area.loc["F277W", "area_arcmin2"], 2.3352655, places=4)
        self.assertAlmostEqual(area.loc["F356W", "area_arcmin2"], 2.3307178, places=4)
        self.assertEqual((int(area.loc["F277W", "peaks_7sigma"]),
                          int(area.loc["F277W", "morphology_pass"])), (148, 105))
        self.assertEqual((int(area.loc["F356W", "peaks_7sigma"]),
                          int(area.loc["F356W", "morphology_pass"])), (103, 80))
        self.assertEqual(len(pd.read_csv(g / "blind_pairs.csv")), 31)
        self.assertEqual(len(pd.read_csv(g / "blind_fixed_pixels.csv")), 28)
        rec = pd.read_csv(g / "blind_recovery.csv")
        for band, k, n in [("F277W", 212, 294), ("F356W", 246, 324)]:
            q = rec[(rec["filter"] == band) & (rec.mag_ab == 27.5)]
            self.assertEqual((int(q.recovered.sum()), int(q.trials.sum())), (k, n))
            self.assertAlmostEqual(q.recovered.sum() / q.trials.sum(),
                                   k / n, places=10)
        dq = pd.read_csv(g / "nuclear_native_checks.csv")
        self.assertEqual(int(dq[dq.site == 185].dq.sum()), 0)

    def test_nexus_joined_selection_and_depth(self):
        n = SRC / "high-z-hosts"
        joint = pd.read_csv(n / "joint_contrasts.csv.gz")
        host = joint[joint["sample"] == "high_z"]
        ctrl = joint[joint["sample"] == "control"]
        self.assertEqual((len(host), host.id.nunique()), (1577, 100))
        self.assertEqual((len(ctrl), int(ctrl.passes.sum())), (5743, 3))
        self.assertEqual(int(host.passes.sum()), 0)
        for band, npos, ntests in [("f200w", 120, 1867), ("f444w", 122, 1936)]:
            c = pd.read_csv(n / band / "contrasts.csv.gz")
            h = c[c["sample"] == "high_z"]
            self.assertEqual((h.id.nunique(), len(h)), (npos, ntests))
        for f, nplace, p90, p50 in [
            ("f200w/recovery.csv", 815, 26.11, 26.70),
            ("f444w/recovery.csv", 839, 26.14, 26.78),
            ("joint_recovery.csv", 686, 26.59, 27.13),
        ]:
            d = pd.read_csv(n / f).sort_values("injected_ab")
            self.assertTrue((d["count"] == nplace).all())
            self.assertTrue(np.allclose(d["mean"], d["sum"] / d["count"]))
            # Interpolate between actual injection bins, never extrapolate.
            for p, expected in [(0.9, p90), (0.5, p50)]:
                i = np.flatnonzero(d["mean"].to_numpy() < p)[0]
                x0, x1 = d.iloc[i - 1], d.iloc[i]
                estimate = x0.injected_ab + (p - x0["mean"]) / (x1["mean"] - x0["mean"]) * (
                    x1.injected_ab - x0.injected_ab)
                self.assertAlmostEqual(estimate, expected, delta=.011)

    def test_cluster_limits_and_archive_receipts(self):
        m = SRC / "macs0416"
        s = json.loads((m / "screen.json").read_text())
        self.assertAlmostEqual(s["area_joint_arcmin2"], .190136, places=6)
        self.assertEqual(s["two_band_detections"], 0)
        inj = pd.read_csv(m / "injections.csv")
        self.assertEqual(int(inj.injected.sum()), 2700)
        q = inj.groupby("ab_each")[["injected", "recovered"]].sum()
        p0 = q.loc[26.5, "recovered"] / q.loc[26.5, "injected"]
        p1 = q.loc[27., "recovered"] / q.loc[27., "injected"]
        mid = 26.5 + .5 * (.5 - p0) / (p1 - p0)
        self.assertAlmostEqual(mid, 26.65, delta=.015)
        lens = json.loads((m / "warhol_lens.json").read_text())
        self.assertAlmostEqual(lens["dist_to_tangential_curve_arcsec"], .410, delta=.001)
        beta = quad(lambda z: 1 / np.sqrt(.3 * (1 + z) ** 3 + .7),
                    .396, .9397)[0] / quad(lambda z: 1 / np.sqrt(.3 * (1 + z) ** 3 + .7),
                                            0, .9397)[0]
        self.assertAlmostEqual(beta, lens["beta"], places=9)
        e = json.loads((SRC / "el-gordo/screen.json").read_text())
        self.assertEqual(e["two_band_detections"], 0)
        # Verify the hashes of deposited inventory bytes, and presence/format
        # of SHA256s for external FITS/catalog bytes (which are not deposited).
        h = SRC / "repeat-catalog-followup/cosmos-022-history"
        entries = pd.read_csv(h / "inventory_checksums.csv")
        self.assertEqual(set(entries.input), {"level2.csv", "level3.csv", "mast_cone.csv"})
        for row in entries.itertuples():
            raw = (h / row.input).read_bytes()
            self.assertEqual((hashlib.sha256(raw).hexdigest(), len(raw)),
                             (row.sha256, row.bytes))
        for file in [
            "repeat-catalog-followup/catalog_inputs.csv",
            "repeat-catalogs/catalog_inputs.csv",
            "goods-s-oct2026/8060_manifest.csv",
            "goods-s-oct2026/12577_manifest.csv",
        ]:
            c = pd.read_csv(SRC / file)
            self.assertTrue(len(c) > 0)
            self.assertTrue(c["url"].str.startswith("https://").all())
            self.assertTrue(c.sha256.str.fullmatch(r"[0-9a-f]{64}").all())
            self.assertTrue((c["bytes"] > 0).all())
        for file in ["macs0416/warhol_lens.json",
                     "repeat-catalog-followup/cosmos-022-validation/muse_receipt.json"]:
            self.assertTrue((SRC / file).stat().st_size > 100)
        receipts = json.loads((m / "warhol_lens.json").read_text())["inputs"]
        self.assertEqual(len(receipts), 2)
        self.assertTrue(all(re.fullmatch(r"[0-9a-f]{64}", x["sha256"]) for x in receipts))


if __name__ == "__main__":
    unittest.main()
