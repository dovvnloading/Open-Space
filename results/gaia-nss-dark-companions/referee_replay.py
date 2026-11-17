"""Independent offline audit of deposited Gaia NSS screens, epoch receipts and fits.

Run from anywhere: python results/gaia-nss-dark-companions/referee_replay.py
Optional: --check-fit /tmp/independent/lamost_candidate_fits.json
The full ESA 167 MB NSS bulk file and 12.5 MB TAP response are *not* deposited.
This tests every deposited screening row and the frozen survey responses; use the
upstream acquisition scripts for a fresh end-to-end source download.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as src:
        return list(csv.DictReader(src))


def load_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_hash(path, expected):
    actual = sha256(path)
    assert actual == expected, (str(path.relative_to(ROOT)), actual, expected)
    return actual


def assert_near(actual, expected, tol=0.00001):
    assert abs(float(actual)-expected) <= tol, (actual, expected, tol)


def verify_rows():
    initial = read_csv(ROOT / "gaia/mass_screen.csv")
    extended = read_csv(ROOT / "extension/mass_screen.csv")
    initial_depth = load_json(ROOT / "gaia/depth.json")
    ext_depth = load_json(ROOT / "extension/depth.json")
    assert len(initial) == initial_depth["thresholds_msun"]["1.4"]["rows"] == 984
    assert len({r["source_id"] for r in initial}) == 983
    assert len(initial) == sum(initial_depth["thresholds_msun"]["1.4"]["by_type"].values())
    assert sum(ext_depth["rows_by_solution"].values()) == len(extended)
    assert all(r["flame_flags"] in ("00", "10") for r in extended)
    seen = set()
    for label, rows in (("initial", initial), ("FLAME", extended)):
        for r in rows:
            key = (r["source_id"], r["nss_solution_type"])
            assert key not in seen or label == "FLAME", (label, key)
            f, m2 = float(r["f_msun"]), float(r["m2_screen_msun"])
            m1 = float(r["m1"]) if r["m1"] else 0.0
            assert f > 0 and m2 > 1.4 and m1 >= 0 and all(map(math.isfinite, [f, m2, m1]))
            reconstructed = m2**3/(m1+m2)**2
            assert abs(reconstructed-f) <= 5e-10*max(1, f), (label, key, reconstructed, f)
            seen.add(key)
    initial_ids = {r["source_id"] for r in initial}
    ext_ids = {r["source_id"] for r in extended}
    high_initial = {r["source_id"] for r in initial if float(r["m2_screen_msun"])>3}
    high_ext = {r["source_id"] for r in extended if float(r["m2_screen_msun"])>3}
    assert len(high_initial)==112
    assert len(ext_ids)==6797
    assert len(high_ext)==291
    assert len(ext_ids-initial_ids)==6750
    assert len(high_ext-high_initial)==271
    assert len(initial_ids|ext_ids)==7733
    assert len(high_initial|high_ext)==383
    assert ext_depth["thresholds_union"] == [7733,383]
    assert ext_depth["new_to_1p4_screen"]==6750
    assert ext_depth["new_to_3_stratum"]==271
    for file,expected in load_json(ROOT/"vetting/seismic_depth.json")["inputs"].items():
        if file in ("extension/mass_screen.csv","gaia/mass_screen.csv"):
            check_hash(ROOT/file,expected["sha256"])
    print("Gaia screen: 984 rows / 983 distinct; 112 above 3; cubic mass equations verified")
    print("FLAME extension: %d rows / 6797 distinct; 6750 novel to cut; 271 new above 3" % len(extended))
    return initial_ids|ext_ids, high_initial|high_ext


def verify_coverage(ids, high_ids):
    meta=load_json(ROOT/"coverage.json")
    with gzip.open(ROOT/"coverage.csv.gz","rt",newline="") as inp:
        coverage=list(csv.DictReader(inp))
    assert len(coverage)==len(ids)==meta["above_1p4"]["sources"]
    assert {r["source_id"] for r in coverage}==ids
    def get_statistics(group):
        names=("lamost","apogee","boss","desi")
        return dict(sources=len(group),
            matched_any_prior_table=sum(int(r["published_study_matches"])>0 for r in group),
            any_clean_survey=sum(any(int(r[k+"_clean"])>0 for k in names) for r in group),
            with_clean_repeats_same_instrument=sum(any(int(r[k+"_max_nights_one_instrument"])>=2 for k in names) for r in group),
            with_three_clean_nights_same_instrument=sum(any(int(r[k+"_max_nights_one_instrument"])>=3 for k in names) for r in group))
    for section,group in (("above_1p4",coverage),("above_3",[r for r in coverage if r["source_id"] in high_ids])):
        computed=get_statistics(group)
        for key,value in computed.items(): assert meta[section][key]==value,(section,key,value,meta[section][key])
        for survey in ("lamost","apogee","boss","desi"):
            assert meta[section][survey]["clean_sources"]==sum(int(r[survey+"_clean"])>0 for r in group)
            assert meta[section][survey]["clean_multinight"]==sum(int(r[survey+"_clean_nights"])>=2 for r in group)
    assert meta["above_1p4"]["any_clean_survey"]==568
    assert meta["above_1p4"]["with_clean_repeats_same_instrument"]==257
    assert meta["above_3"]["any_clean_survey"]==24
    assert meta["above_3"]["with_clean_repeats_same_instrument"]==12
    assert meta["above_3"]["sources"]-24==359
    refs=set()
    for folder in ("novelty","extension/novelty"):
        for r in read_csv(ROOT/folder/"prior_matches.csv"):
            assert r["source_id"] in ids
            refs.add(r["source_id"])
    assert len(refs)==1627
    assert sum(int(r["published_study_matches"])>0 for r in coverage)==len(refs)
    print("Coverage: 7733 sources / 568 clean / 257 repeated; 383 high / 24 clean / 12 repeated")
    print("Exact-ID prior tables: 1627 selected sources; absence is NOT a novelty claim")


def verify_receipts():
    stats={}
    for folder in ("lamost","extension/lamost"):
        root=ROOT/folder/"raw"
        receipts=load_json(root/"receipts.json")
        packed=root/"responses.zip"
        z=ZipFile(packed) if packed.exists() else None
        try:
            total_rows=0
            for r in receipts:
                filename=r["file"]
                data=(root/filename).read_bytes() if (root/filename).exists() else z.read(filename)
                raw=gzip.decompress(data)
                h=hashlib.sha256(raw).hexdigest()
                assert h==r["sha256_csv"],(folder,filename,h,r["sha256_csv"])
                parsed=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
                assert len(parsed)==r["rows"],(folder,filename,len(parsed),r["rows"])
                requested=set(r["requested_ids"])
                assert all(p["gaia_source_id"] in requested for p in parsed)
                total_rows+=len(parsed)
            stats[folder]=dict(receipts=len(receipts),returned_rows=total_rows)
        finally:
            if z:z.close()
    assert stats["lamost"]["returned_rows"]==load_json(ROOT/"lamost/depth.json")["returned_rows"]
    assert stats["extension/lamost"]["returned_rows"]==load_json(ROOT/"extension/lamost/depth.json")["returned_rows"]
    print("LAMOST original/extension compressed server replies verified:",stats)
    return stats


def verify_dynamics():
    candidate=load_json(ROOT/"orbit/lamost_candidate_fits.json")
    source_hashes=candidate["sha256"]
    for name,path in {"lamost_candidate_nss.csv":ROOT/"orbit/lamost_candidate_nss.csv",
                      "exposures.csv":ROOT/"lamost/exposures.csv",
                      "mass_screen.csv":ROOT/"gaia/mass_screen.csv"}.items():
        check_hash(path,source_hashes[name])
    exposures=[r for r in read_csv(ROOT/"lamost/exposures.csv") if r["source_id"]==candidate["source_id"]]
    assert len(exposures)==candidate["exposures"]==100
    assert len({r["obsdate"] for r in exposures})==18
    primary=candidate["models"]["NSSx1_floor1_M11.950"]
    assert_near(primary["period_d"][1],1007.3,.1)
    assert_near(primary["m2_msun"][1],1.77,.003)
    assert_near(primary["m2_msun"][0],1.71,.005)
    assert_near(primary["m2_msun"][2],1.84,.005)
    assert_near(primary["rv_rms_kms"],.69,.005)
    assert_near(primary["f_rv_msun"][1],.384,.002)
    assert_near(candidate["holdout_fit"]["holdout"]["rms_kms"],1.20,.01)
    residuals=[r for r in read_csv(ROOT/"orbit/lamost_candidate_residuals.csv")
               if r["fit"]=="NSSx1_floor1_M11.950"]
    assert len(residuals)==18
    rms=math.sqrt(sum(float(r["observed_minus_model_kms"])**2 for r in residuals)/len(residuals))
    assert_near(rms,primary["rv_rms_kms"],.000005)
    giant=load_json(ROOT/"orbit/giant_sb1_fits.json")
    for name,path in {"flame_giant_nss.csv":ROOT/"orbit/flame_giant_nss.csv",
                      "exposures.csv":ROOT/"extension/lamost/exposures.csv",
                      "mass_screen.csv":ROOT/"extension/mass_screen.csv"}.items():
        check_hash(path,giant["input_sha256"][name])
    g=giant["results"]["nominal"]
    assert giant["n_nights"]==7 and giant["n_exposures"]==30
    assert_near(g["mass_function_msun"][1],.640,.002)
    assert_near(g["critical_primary_msun_for_min_m2_3"][1],3.49,.02)
    assert_near(g["mass_models"]["catalog"]["m2_min_q16_q50_q84"][1],3.60,.02)
    assert_near(g["mass_models"]["30pct_lower_primary"]["m2_min_q16_q50_q84"][1],2.99,.02)
    desi=load_json(ROOT/"orbit/desi_lpv_fits.json")
    assert_near(desi["span_days"],45.87,.02)
    assert_near(desi["models"]["nominal"]["f_rv_msun"][1],.561,.008)
    controls=load_json(ROOT/"orbit/control_fits.json")
    for name,mass in {"BH1":9.63,"BH2":8.93,"NS1":1.90,"BH3":33.12}.items():
        assert_near(controls[name]["m2_msun"][1],mass,.02)
    print("HD 264291: fit-file/input hashes OK; 18 residual nights, mass 1.77, RMS 0.69, holdout 1.20")
    print("Giant: 470-day orbit, f 0.640, M2 floor 3.60 vs 2.99 with lower M1")
    print("DESI LPV and four known compact controls: deposited fits and values checked")


def check_fresh_fit(path):
    actual=load_json(path)
    saved=load_json(ROOT/"orbit/lamost_candidate_fits.json")
    assert actual["sha256"]==saved["sha256"]
    assert actual["exposures"]==saved["exposures"]==100
    assert actual["source_id"]==saved["source_id"]
    for label in saved["models"]:
        a,b=actual["models"][label],saved["models"][label]
        assert len(a["m2_msun"])==len(b["m2_msun"])==3
        for metric,tolerance in (("period_d",.45),("m2_msun",.018),("f_rv_msun",.009),
                                 ("inclination_deg",.65)):
            for x,y in zip(a[metric],b[metric]):
                assert_near(x,y,tolerance)
        assert_near(a["rv_rms_kms"],b["rv_rms_kms"],.002)
        assert_near(a["rv_chi2"],b["rv_chi2"],.02)
    ah,bh=actual["holdout_fit"]["holdout"],saved["holdout_fit"]["holdout"]
    assert ah["n"]==bh["n"]==7
    assert_near(ah["rms_kms"],bh["rms_kms"],.12)
    assert_near(ah["predictive_tail_fraction"],bh["predictive_tail_fraction"],.08)
    print("FRESH HD 264291 FIT RECOVERED: mass",actual["models"]["NSSx1_floor1_M11.950"]["m2_msun"],
          "period",actual["models"]["NSSx1_floor1_M11.950"]["period_d"],
          "holdout RMS",ah["rms_kms"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-fit",type=Path,help="Compare a fresh fit in a separate output folder")
    a=parser.parse_args()
    if a.check_fit:
        check_fresh_fit(a.check_fit)
    else:
        ids,high=verify_rows()
        verify_coverage(ids,high)
        verify_receipts()
        verify_dynamics()
        print("Offline referee check PASS. Raw ESA bulk/TAP source acquisition is a separate upstream replay.")


if __name__=="__main__":
    main()
