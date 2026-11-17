"""Independent ledger-and-receipt referee tests for JWST five-field PR #146.

Read committed outputs, not README prose or other test modules. This is an
offline cross-check; FITS-level reprocessing additionally needs the hashed
multi-GB MAST files specified in each field's input manifest.
"""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parents[1] / "results/jwst-five-field-followup"
L = HERE / "lawd66"
N = HERE / "nircam-hd53143-hd38206/outputs"
D = HERE / "hd21997-hd202628/measurements"
M = HERE / "hd53143-miri"

def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

class IndependentReferee(unittest.TestCase):
    def test_lawd_registration_and_unresolved_peak(self):
        epochs = sorted((x for x in read(L / "epochs.csv") if x["filter"] == "F164N"), key=lambda x: float(x["mjd"]))
        self.assertEqual(len(epochs), 8)
        a, z = epochs[0], epochs[-1]
        cosdec = math.cos(math.radians(float(a["dec"])))
        raw_e = 3600 * (float(z["ra"]) - float(a["ra"])) * cosdec
        raw_n = 3600 * (float(z["dec"]) - float(a["dec"]))
        registered_e = raw_e + float(z["east_tie_shift"]) - float(a["east_tie_shift"])
        registered_n = raw_n + float(z["north_tie_shift"]) - float(a["north_tie_shift"])
        self.assertAlmostEqual(math.hypot(registered_e, registered_n), 1.207, places=3)
        self.assertAlmostEqual(float(z["mjd"]) - float(a["mjd"]), 1140.92, places=2)
        peaks = read(L / "peaks.csv")
        groups = read(L / "sources.csv")
        self.assertEqual((len(peaks), len(groups)), (93, 10))
        self.assertEqual(sum(int(p["fixed_wins"]) > 0 for p in peaks), 92)
        self.assertEqual(sum(int(g["fixed_wins"]) > 0 for g in groups), 9)
        self.assertEqual(sum(int(p["host_wins"]) > 0 for p in peaks), 0)
        lone = [p for p in peaks if not int(p["fixed_wins"])]
        self.assertEqual(len(lone), 1)
        p = lone[0]
        self.assertEqual(p["band"], "F164N")
        self.assertAlmostEqual(float(p["r_arcsec"]), 3.274, places=3)
        self.assertAlmostEqual(float(p["east_arcsec"]), 0.344, places=3)
        self.assertAlmostEqual(float(p["north_arcsec"]), 3.256, places=3)
        self.assertAlmostEqual(float(p["score"]), 6.10, places=2)
        d = read(L / "dither_photometry.csv")
        self.assertEqual(len(d), 4)
        self.assertTrue(all(float(x["forced_flux_uJy"]) > 0 and int(x["aperture_DQ_bitwise_or"]) == 0 for x in d))

    def test_hd53143_background_astrometry_color_and_miri_null(self):
        models = read(N / "NICI_track.csv")
        self.assertEqual(len(models), 2)
        for m in models:
            x, y = float(m["measured_east"]), float(m["measured_north"])
            self.assertAlmostEqual(math.hypot(x-float(m["pred_bg_east"]),y-float(m["pred_bg_north"])),float(m["residual_bg"]))
            self.assertAlmostEqual(float(m["residual_bg"]), 0.149, places=3)
            self.assertGreater(float(m["residual_comoving"]), 4.67)
        colours = read(N / "known_source_color.csv")
        self.assertEqual(len(colours), 2)
        self.assertTrue(all(abs(float(x["color_AB"]) + 1.17) < 0.01 for x in colours))
        fits = list((M/"measurements").glob("*/peaks.csv"))
        self.assertEqual(len(fits), 3)
        self.assertTrue(all(len(read(f)) == 0 for f in fits))

    def test_hd38206_and_other_disk_fields(self):
        jan, octo = read(N/"hd38206_o141_peaks.csv"), read(N/"hd38206_o150_peaks.csv")
        self.assertEqual((len(jan),len(octo)), (5,3))
        repeat = {(p["x0"], p["y0"]) for p in jan} & {(p["x0"], p["y0"]) for p in octo}
        self.assertEqual(repeat, {("98","215"),("166","235")})
        tracks = [p for p in read(N/"paired_rolls.csv") if p["origin"].startswith("hd38206") and p["hypothesis"] == "comoving"]
        self.assertEqual(len(tracks), 8)
        self.assertLess(max(float(p["max_3x3_snr"]) for p in tracks), 2.36)
        self.assertEqual(len(read(D/"HD21997_peaks.csv")), 0)
        self.assertEqual(len(read(D/"HD202628_peaks.csv")), 2)
        refs = sorted(round(float(p["matched_snr"]),2) for p in read(D/"HD202628_reference_test.csv"))
        self.assertEqual(refs, [-0.50,0.31,5.44,8.16])
        old = read(N/"gemini_locus_summary.csv")
        self.assertEqual(len(old), 16)
        self.assertTrue(all(int(r["measured"]) == 20 for r in old))

    def test_hd202628_f1140c_two_reference_sensitivity(self):
        by_method = defaultdict(lambda: defaultdict(dict))
        for r in read(D/"F1140C_recovery_trials.csv"):
            if r["observed_PA"] == "0" and r["eligible"] == "1":
                by_method[(r["method"], float(r["flux_uJy"]))][r["roll"]][r["trial"]] = int(r["recovered"])
        minima = {}
        for meth in ("native_RDI_HD201865", "native_RDI_HD22333"):
            passing=[]
            for (method,level), rolls in by_method.items():
                if method != meth: continue
                common=rolls["007"].keys() & rolls["008"].keys()
                paired=[rolls["007"][k] and rolls["008"][k] for k in common]
                self.assertEqual(len(paired),28)
                if level == 10: self.assertEqual(sum(paired),9 if "201865" in meth else 0)
                if sum(paired)/len(paired) >= .95: passing.append(level)
            minima[meth]=min(passing)
        self.assertEqual(minima, {"native_RDI_HD201865":15., "native_RDI_HD22333":20.})
        # Different native references disagree at a flux smaller than both limits.
        rows = read(D/"MIRI_roll_tests.csv")
        self.assertTrue(any(r["method"] == "pipeline" and abs(float(r["fixed_snr"])-6.1202848455)<.001 for r in rows))
        self.assertTrue(any(r["method"].startswith("native") and abs(float(r["fixed_snr"])+.77)<.10 for r in rows))

    def test_all_57_empirical_contrast_points_against_individual_trials(self):
        summary=read(HERE/"contrast_summary.csv")
        trials=[]
        for visit,field,series in [
            ("hd53143_o001","HD 53143","Roll 001"),("hd53143_o008","HD 53143","Roll 008"),
            ("hd38206_o141","HD 38206","Jan 2025"),("hd38206_o150","HD 38206","Oct 2025")]:
            trials.append((read(N/f"{visit}_injections.csv"),field,"F444W",series,"contrast",.95,True))
        for suffix,series in [("21997","Feb 2024"),("202628","Jun 2024")]:
            trials.append((read(D/f"HD{suffix}_injections.csv"),f"HD {suffix}","F444W",series,"contrast",.95,True))
        for visit,series in [("o003","Roll 003"),("o004","Roll 004"),("c1005","Combined")]:
            trials.append((read(M/f"measurements/jw08142-{visit}/injection_trials.csv"),"HD 53143","F2550W",series,"contrast",.95,False))
        for band in ("F164N","F277W","F323N"):
            trials.append(([r for r in read(L/"joint_injections.csv") if r["band"]==band],"LAWD 66",band,band,"ratio",.90,False))
        found=0
        total=0
        for ledger,field,band,series,key,goal,sustained in trials:
            total+=len(ledger)
            grouped=defaultdict(lambda: defaultdict(list))
            for item in ledger:
                grouped[float(item["radius_arcsec"])][float(item[key])].append(int(item["recovered"]))
            for radius,levels in grouped.items():
                pct={v:sum(a)/len(a) for v,a in levels.items()}
                good=sorted(v for v in levels if pct[v]>=goal and
                    (not sustained or all(pct[w]>=goal for w in levels if w>=v)))
                self.assertTrue(good)
                expected=next(s for s in summary if
                    (s["field"],s["band"],s["series"],float(s["radius_arcsec"]))==(field,band,series,radius))
                limit=good[0]
                self.assertTrue(math.isclose(limit,float(expected["limiting_tested_flux_ratio"]),rel_tol=1e-12))
                self.assertEqual((sum(levels[limit]),len(levels[limit])),
                                 (int(expected["recovered_at_limit"]),int(expected["tried_at_limit"])))
                found+=1
        self.assertEqual((found,len(summary)),(57,57))
        self.assertEqual(total+len(read(L/"injections.csv")),26809)

    def test_archival_receipts_and_hst_controls(self):
        manifest=json.loads((L/"inputs.json").read_text())
        for row in manifest["mast_queries"]+manifest["receipts"]:
            self.assertEqual(hashlib.sha256((L/row["file"]).read_bytes()).hexdigest(),row["sha256"])
        query=json.loads((M/"receipts/query.json").read_text())
        self.assertEqual(hashlib.sha256((M/"receipts/caom.json").read_bytes()).hexdigest(),query["sha256"])
        a=HERE/"archives"
        for q in json.loads((a/"query_receipts.json").read_text()):
            raw=gzip.decompress((a/"receipts"/(q["target"]+".csv.gz")).read_bytes())
            self.assertEqual(hashlib.sha256(raw).hexdigest(),q["uncompressed_sha256"])
        self.assertEqual(Counter(r["field"] for r in read(a/"hst_tracks.csv")),
                         {"HD53143":504,"HD38206":32,"HD21997":12,"HD202628":168})

if __name__ == "__main__":
    unittest.main()
