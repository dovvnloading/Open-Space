"""Replay the five-field dispositions from committed per-source and trial ledgers."""

import csv
import gzip
import hashlib
import json
import math
import unittest
from collections import Counter, defaultdict
from pathlib import Path

RESULT = Path(__file__).resolve().parents[1] / "results/jwst-five-field-followup"
NIR = RESULT / "nircam-hd53143-hd38206/outputs"
DISKS = RESULT / "hd21997-hd202628/measurements"
LAWD = RESULT / "lawd66"
MIRI = RESULT / "hd53143-miri"


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


class FiveFieldFollowup(unittest.TestCase):
    def test_lawd_positions_and_tentative_peak(self):
        groups = rows(LAWD / "sources.csv")
        peaks = rows(LAWD / "peaks.csv")
        self.assertEqual((len(groups), len(peaks)), (10, 93))
        ids = [peak for group in groups for peak in group["peak_ids"].split()]
        self.assertCountEqual(ids, [row["id"] for row in peaks])
        fixed = {row["id"] for row in peaks if int(row["fixed_wins"]) > 0}
        self.assertEqual(len(fixed), 92)
        self.assertEqual(
            sum(bool(fixed.intersection(g["peak_ids"].split())) for g in groups), 9
        )
        self.assertTrue(all(int(row["host_wins"]) == 0 for row in peaks))
        dithers = rows(LAWD / "dither_photometry.csv")
        self.assertEqual(len(dithers), 4)
        self.assertTrue(
            all(
                float(r["forced_flux_uJy"]) > 0
                and int(r["aperture_DQ_bitwise_or"]) == 0
                for r in dithers
            )
        )
        native = rows(LAWD / "unmatched_forced.csv")
        for fit in rows(LAWD / "unmatched_constant_flux.csv"):
            subset = [
                r
                for r in native
                if r["band"] == "F164N"
                and r["prediction"] == fit["prediction"]
                and r["footprint_valid"] == "1"
                and (fit["trigger_epoch_included"] == "1" or r["date"] != "2025-04-21")
            ]
            x = [float(r["flux_uJy"]) for r in subset]
            w = [float(r["err_uJy"]) ** -2 for r in subset]
            avg = sum(v * k for v, k in zip(x, w)) / sum(w)
            self.assertAlmostEqual(avg, float(fit["constant_flux_uJy"]))
            self.assertAlmostEqual(
                sum(k * (v - avg) ** 2 for v, k in zip(x, w)), float(fit["raw_chi2"])
            )

    def test_disk_field_tracks_and_reference_substitutions(self):
        for row in rows(NIR / "NICI_track.csv"):
            com_e = float(row["published_rho"]) * math.sin(
                math.radians(float(row["published_pa"]))
            )
            com_n = float(row["published_rho"]) * math.cos(
                math.radians(float(row["published_pa"]))
            )
            e, n = float(row["measured_east"]), float(row["measured_north"])
            self.assertAlmostEqual(
                math.hypot(e - com_e, n - com_n), float(row["residual_comoving"])
            )
            self.assertAlmostEqual(
                math.hypot(
                    e - float(row["pred_bg_east"]), n - float(row["pred_bg_north"])
                ),
                float(row["residual_bg"]),
            )
            self.assertLess(float(row["residual_bg"]), 0.15)
            self.assertGreater(float(row["residual_comoving"]), 4.67)
        tracks = [
            r
            for r in rows(NIR / "paired_rolls.csv")
            if r["origin"].startswith("hd38206")
        ]
        fixed = [r for r in tracks if r["hypothesis"] == "comoving"]
        self.assertEqual(len(fixed), 8)
        self.assertLess(max(float(r["max_3x3_snr"]) for r in fixed), 2.36)
        self.assertEqual(len(rows(DISKS / "HD21997_peaks.csv")), 0)
        self.assertEqual(len(rows(DISKS / "HD202628_peaks.csv")), 2)
        refs = rows(DISKS / "HD202628_reference_test.csv")
        self.assertEqual(
            sorted(round(float(r["matched_snr"]), 2) for r in refs),
            [-0.50, 0.31, 5.44, 8.16],
        )
        self.assertEqual(len(rows(DISKS / "all_model_comparisons.csv")), 17)
        files = list((MIRI / "measurements").glob("*/peaks.csv"))
        self.assertEqual(len(files), 3)
        self.assertTrue(all(len(rows(f)) == 0 for f in files))

    def test_recovery_limits_from_individual_trials(self):
        expected = rows(RESULT / "contrast_summary.csv")
        datasets = [
            (
                NIR / f"{visit}_injections.csv",
                field,
                "F444W",
                series,
                "contrast",
                0.95,
                True,
            )
            for visit, field, series in [
                ("hd53143_o001", "HD 53143", "Roll 001"),
                ("hd53143_o008", "HD 53143", "Roll 008"),
                ("hd38206_o141", "HD 38206", "Jan 2025"),
                ("hd38206_o150", "HD 38206", "Oct 2025"),
            ]
        ]
        datasets += [
            (
                DISKS / f"HD{hd}_injections.csv",
                f"HD {hd}",
                "F444W",
                series,
                "contrast",
                0.95,
                True,
            )
            for hd, series in [("21997", "Feb 2024"), ("202628", "Jun 2024")]
        ]
        datasets += [
            (
                MIRI / f"measurements/jw08142-{visit}/injection_trials.csv",
                "HD 53143",
                "F2550W",
                series,
                "contrast",
                0.95,
                False,
            )
            for visit, series in [
                ("o003", "Roll 003"),
                ("o004", "Roll 004"),
                ("c1005", "Combined"),
            ]
        ]
        datasets += [
            (LAWD / "joint_injections.csv", "LAWD 66", band, band, "ratio", 0.90, False)
            for band in ("F164N", "F277W", "F323N")
        ]
        total = 0
        checked = 0
        for path, field, band, series, key, goal, sustained in datasets:
            trials = [r for r in rows(path) if r.get("band", band) == band]
            total += len(trials)
            by_radius = defaultdict(lambda: defaultdict(list))
            for r in trials:
                by_radius[float(r["radius_arcsec"])][float(r[key])].append(
                    int(r["recovered"])
                )
            for radius, levels in by_radius.items():
                rec = {k: sum(v) / len(v) for k, v in levels.items()}
                good = [
                    k
                    for k in levels
                    if rec[k] >= goal
                    and (
                        not sustained
                        or all(
                            value >= goal for test, value in rec.items() if test >= k
                        )
                    )
                ]
                reported = next(
                    r
                    for r in expected
                    if (r["field"], r["band"], r["series"], float(r["radius_arcsec"]))
                    == (field, band, series, radius)
                )
                self.assertTrue(
                    math.isclose(
                        min(good),
                        float(reported["limiting_tested_flux_ratio"]),
                        rel_tol=1e-12,
                    )
                )
                self.assertEqual(
                    sum(levels[min(good)]), int(reported["recovered_at_limit"])
                )
                self.assertEqual(
                    len(levels[min(good)]), int(reported["tried_at_limit"])
                )
                checked += 1
        self.assertEqual(checked, len(expected))
        self.assertEqual(total + len(rows(LAWD / "injections.csv")), 26809)

    def test_f1140c_joint_recovery_at_unclassified_peak(self):
        by_level = defaultdict(lambda: defaultdict(dict))
        for row in rows(DISKS / "F1140C_recovery_trials.csv"):
            if row["observed_PA"] == "0" and row["eligible"] == "1":
                key = (row["method"], float(row["flux_uJy"]))
                by_level[key][row["roll"]][row["trial"]] = int(row["recovered"])
        limits = defaultdict(list)
        for row in rows(DISKS / "F1140C_recovery.csv"):
            key = (row["method"], float(row["flux_uJy"]))
            samples = by_level[key]
            if row["roll"] == "both":
                common = samples["007"].keys() & samples["008"].keys()
                trials = [samples["007"][n] and samples["008"][n] for n in common]
            else:
                trials = list(samples[row["roll"]].values())
            self.assertEqual(
                (sum(trials), len(trials)),
                (int(row["recovered"]), int(row["eligible"])),
            )
            if row["roll"] == "both" and sum(trials) >= 0.95 * len(trials):
                limits[row["method"]].append(key[1])
        self.assertEqual(min(limits["native_RDI_HD201865"]), 15)
        self.assertEqual(min(limits["native_RDI_HD22333"]), 20)

    def test_ground_loci(self):
        grouped = rows(NIR / "gemini_locus_summary.csv")
        self.assertEqual(len(grouped), 16)
        self.assertEqual({int(r["measured"]) for r in grouped}, {20})
        self.assertLess(max(abs(float(r["median_raw_score"])) for r in grouped), 2.26)

    def test_archival_receipts(self):
        manifest = json.loads((LAWD / "inputs.json").read_text())
        for r in manifest["mast_queries"] + manifest["receipts"]:
            self.assertEqual(
                hashlib.sha256((LAWD / r["file"]).read_bytes()).hexdigest(), r["sha256"]
            )
        rr = json.loads((MIRI / "receipts/query.json").read_text())
        self.assertEqual(
            hashlib.sha256((MIRI / "receipts/caom.json").read_bytes()).hexdigest(),
            rr["sha256"],
        )
        archive = RESULT / "archives"
        for q in json.loads((archive / "query_receipts.json").read_text()):
            raw = gzip.decompress(
                (archive / "receipts" / (q["target"] + ".csv.gz")).read_bytes()
            )
            self.assertEqual(hashlib.sha256(raw).hexdigest(), q["uncompressed_sha256"])
        tracks = rows(archive / "hst_tracks.csv")
        self.assertEqual(
            Counter(r["field"] for r in tracks),
            {"HD53143": 504, "HD38206": 32, "HD21997": 12, "HD202628": 168},
        )
