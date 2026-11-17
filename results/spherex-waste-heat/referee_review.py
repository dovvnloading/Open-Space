"""Independent referee replay. Run in an isolated GitHub Actions checkout.

offline: rerun reductions from the deposited bytes, check archive receipts and
numerical claims; fresh: independently read the DR1 FITS tiles and compare
the measured spectra to the frozen tables. No branch mutations are published.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile

import numpy as np

STUDY=Path(__file__).resolve().parent
REPO=STUDY.parents[1]
SUMMARY={}

def rows(path):
    with Path(path).open(newline="") as f: return list(csv.DictReader(f))

def value(x):
    try: return float(x)
    except (TypeError,ValueError): return None

def approx(got,expected,atol=0.006):
    if got is None or not math.isclose(float(got),expected,abs_tol=atol):
        raise AssertionError(f"expected {expected}, received {got}")

def cmd(*args):
    print("+ "+" ".join(map(str,args)),flush=True)
    subprocess.run([str(x) for x in args],cwd=REPO,check=True)

def snapshot(relatives):
    return {r:(REPO/r).read_bytes() for r in relatives}

def compare(saved):
    for rel,previous in saved.items():
        current=(REPO/rel).read_bytes()
        if current!=previous:
            # Decimal serialization can differ between Python/Numpy patch versions.
            if rel.endswith(".json"):
                assert json.loads(current)==json.loads(previous),rel
            elif rel.endswith(".csv"):
                def table(raw):
                    return list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
                a,b=table(previous),table(current)
                assert len(a)==len(b),f"{rel}: {len(a)} != {len(b)}"
                for i,(x,y) in enumerate(zip(a,b)):
                    assert x.keys()==y.keys(),rel
                    for key in x:
                        if x[key]==y[key]:continue
                        vx,vy=value(x[key]),value(y[key])
                        assert vx is not None and vy is not None and math.isclose(
                            vx,vy,abs_tol=1e-5,rel_tol=2e-5
                        ),f"{rel} row {i} {key}: {vx} -> {vy}"
            else: raise AssertionError(f"Nonreproduced {rel}")
    return len(saved)

def claim_audit():
    import collections
    blind=json.loads((STUDY/"blind/counts.json").read_text())
    dr1=json.loads((STUDY/"blind/dr1_counts.json").read_text())
    fgk=json.loads((STUDY/"fgk_spectra/counts.json").read_text())
    vet=json.loads((STUDY/"blind/classifications/counts.json").read_text())
    assert blind["field"]==596 and blind["screened_either"]==404
    assert blind["screened_w3"]==404 and blind["screened_w4"]==324
    assert not blind["photometric_candidates"] and blind["injection_controls"]==334
    assert blind["first_failure_counts"]=={
       "Gaia astrometry or CMD":115,"no AllWISE match":19,
       "nonunique match":5,"W1 or morphology quality":53,"eligible":404}
    assert vet["complete"] and vet["wise"]["retain_after_secure_and_ambiguous"]==397
    assert vet["spherex"]["retain_after_secure_and_ambiguous"]==54
    assert dr1["strict"]["attempted"]==58 and dr1["strict"]["valid_target_photometry"]==56
    assert dr1["strict"]["tested"]==55 and dr1["strict"]["injection_controls"]==43
    assert dr1["outlier_included"]["tested"]==55 and not dr1["triggers_both_masks"]
    assert fgk["published_target_count"]==53
    assert fgk["classification_counts"]=={
        "redder_than_both_in_both_masks":5,
        "below_threshold":13,
        "insufficient_controls_or_coverage":35}
    fs=rows(STUDY/"fgk_spectra/all_53.csv")
    assert len(fs)==53 and len(set(r["source_id"] for r in fs))==53
    ledger=rows(STUDY/"literature/fgk_dispositions.csv")
    assert len(ledger)==53 and {r["gaia_dr3_id"] for r in ledger}=={r["source_id"] for r in fs}
    assert len([r for r in fs if float(r["distance_inv_parallax_pc"])<=300])==2
    for x in ledger:
        assert x["primary_references"].startswith("https://doi.org/10.3847/1538-3881/ad6b90")
        assert x["spherex_dr1_disposition"] and x["prior_identification_or_limit"]
    literature=(STUDY/"literature/README.md").read_text()
    for letter in "ABCDEFGHIJ":
        assert re.search(r"^\| "+letter+r" \|",literature,re.M),letter
    wise=rows(STUDY/"blind/injections.csv")
    dr_inj=rows(STUDY/"blind/dr1_injections.csv")
    def recovered(table,t,g,mask=None):
        matches=[r for r in table if int(r["dust_K"])==t and
                 math.isclose(float(r["L_excess_over_L_star"]),g) and
                 (mask is None or r["extraction"]==mask)]
        assert len(matches)==1
        k="recovered" if mask else "recovered_either"
        return int(matches[0][k]),int(matches[0].get("injected",matches[0].get("controls")))
    assert [recovered(wise,t,.01) for t in (100,300,600)]==[(84,334),(202,334),(157,334)]
    assert [recovered(wise,t,.1) for t in (100,300,600)]==[(232,334),(334,334),(330,334)]
    assert [recovered(dr_inj,t,.01,"strict") for t in (100,300,600)]==[(0,43),(0,43),(25,43)]
    assert [recovered(dr_inj,t,.1,"strict") for t in (100,300,600)]==[(0,43),(22,43),(43,43)]
    assert round(100*25/43)==58
    screen=rows(STUDY/"blind/screen.csv")
    assert sum(r["status"]=="screened" for r in screen)==404
    assert sum(int(r["selected3"] or 0) or int(r["selected4"] or 0) for r in screen)==0
    assert sum(value(r["snr3"]) is not None and value(r["snr3"])>=5 for r in screen)==205
    assert sum(value(r["snr4"]) is not None and value(r["snr4"])>=5 for r in screen)==0
    w3={r["label"]:r for r in rows(STUDY/"imaging/centroids.csv") if r["band"]=="3"}
    for star,shift in [("B",3.061),("C",3.284),("G",5.539)]:
        approx(w3[star]["separation_as"],shift,1e-5)
    warm={r["label"]:r for r in rows(STUDY/"fgk_spectra/scenes/centroids.csv") if r["band"]=="3"}
    approx(warm["FGK31"]["gaia_offset_arcsec"],.665,1e-5)
    visit=rows(STUDY/"wise/variability.csv")
    assert sum(int(r["raw_neighbors"]) for r in visit if r["band"]=="W1")==173
    assert [int(r["visits_ge3"]) for r in visit if r["label"]=="A"]==[2,2]
    assert all(r["label"]=="A" or int(r["visits_ge3"])==1 for r in visit)
    for r in visit:
        if r["label"]=="A":approx(r["signed_later_minus_earlier_mag"],.007 if r["band"]=="W1" else .0055,.0001)
    stats=json.loads((STUDY/"literature/fgk_statistics.json").read_text())
    assert stats["clean_catwise"]["n"]==53 and stats["clean_catwise"]["n_both_brighter_0p10"]==53
    assert stats["neo_stars_with_visit"]==[52,46] and stats["neo_stars_with_two_visits"]==[4,3]
    curves=rows(STUDY/"warm_fgk/color_grid.csv")
    for sid,pred in [("3017430205815280384",5.07),("6359523376149191936",1.88)]:
        central=[r for r in curves if r["source_id"]==sid and r["catalog"]=="allwise"
                 and float(r["stellar_K"])==5500 and float(r["dW1"])==0 and float(r["dW2"])==0]
        assert len(central)==1
        approx(central[0]["W3_observed_over_T12"],pred,.007)
    skymapper=json.loads((STUDY/"warm_fgk/skymapper/summary.json").read_text())
    assert skymapper["complete"] and all(x["nights_at_least_0p2mag_and_3sigma"]==0 for x in skymapper["series"])
    report=(STUDY/"README.md").read_text()
    for quote in ["0.665", "708 pc", "529 pc", "0/404", "0/55", "0/397", "0/54"]:
        assert quote in report,quote
    assert ("[REMARKABLE]" not in report)
    root=(REPO/"README.md").read_text()
    block=root.split("## All results by topic")[1].split("## How these results were made")[0]
    entries=[r for r in block.splitlines() if r.startswith("- ") and "*Null result*" in r]
    assert len(entries)==18 and "192 entries" in root
    assert sum("results/spherex-waste-heat/README.md" in r for r in entries)==1
    SUMMARY["numbers"]={"WISE":{"selected":0,"eligible":404,"after_blend_veto":397},
      "SPHEREx":{"selected":0,"tested":55,"after_blend_veto":54},
      "WISE_injection_1pct":[84,202,157],
      "SPHEREx_injection_1pct":[0,0,25],"SPHEREx_controls":43,
      "FGK_classifications":[5,13,35],"per_source_FGK_dispositions":53}
    return True

def offline():
    claim_audit()
    paths=[
      "blind/screen.csv","blind/comparison_ids.csv","blind/counts.json",
      "blind/injections.csv","blind/injection_controls.csv",
      "blind/dr1_screen.csv","blind/dr1_counts.json","blind/dr1_injections.csv",
      "literature/fgk_catwise_comparison.csv","literature/fgk_neo_visits.csv",
      "literature/fgk_neo_coverage.csv","literature/fgk_statistics.json",
      "fgk_spectra/all_53.csv","fgk_spectra/counts.json",
      "wise/catalog_photometry.csv","wise/visits.csv","wise/variability.csv",
      "warm_fgk/color_grid.csv","warm_fgk/color_input_hashes.json"
    ]
    original=snapshot(["results/spherex-waste-heat/"+p for p in paths])
    with tempfile.TemporaryDirectory() as tmp:
        archive=Path(tmp)/"wise_archive"
        archive.mkdir()
        with zipfile.ZipFile(STUDY/"wise/archive_inputs.zip") as z: z.extractall(archive)
        output=Path(tmp)/"wise_derived"
        cmd(sys.executable,STUDY/"wise/measure.py","--archive",archive,"--output",output,"--allow-legacy-status")
        for name in ("catalog_photometry.csv","visits.csv","variability.csv"):
            old=(STUDY/"wise"/name).read_bytes()
            new=(output/name).read_bytes()
            assert old==new,f"WISE replay differs: {name}"
        receipt=json.loads((output/"replay-provenance.json").read_text())
        SUMMARY["wise_archive"]={"responses":40,"unknown_status":len(receipt["unknown_status_receipts"])}
    cmd(sys.executable,STUDY/"blind/screen.py")
    cmd(sys.executable,STUDY/"blind/dr1_screen.py")
    cmd(sys.executable,STUDY/"literature/fgk_measure.py")
    cmd(sys.executable,STUDY/"fgk_spectra/run.py","summarize","--prior-controls",STUDY/"spectra/comparison_wise.csv")
    cmd("node",STUDY/"warm_fgk/color_check.mjs")
    SUMMARY["replayed_derived_tables"]=compare(original)
    claim_audit()
    SUMMARY["status"]="passed"
    (STUDY/"referee_offline.json").write_text(json.dumps(SUMMARY,indent=2,sort_keys=True)+"\n")
    print(json.dumps(SUMMARY,indent=2),flush=True)

def fresh():
    import warnings
    import sys as _sys
    spectra=STUDY/"spectra"
    first=["channel_spectra.csv","feature_tests.csv","continuum.csv","coverage.csv"]
    old=snapshot(["results/spherex-waste-heat/spectra/"+p for p in first])
    receipt_dir=spectra/"cube_receipts"
    receipts={p.stem:json.loads(p.read_text()) for p in receipt_dir.glob("*.json")}
    assert len(receipts)==16, f"Expected 16 first-sample cubes, got {len(receipts)}"
    with tempfile.TemporaryDirectory() as tmp:
        cache=Path(tmp)/"stamps"
        for group in ("positions.csv","positions_fgk.csv","positions_pilot.csv"):
            cmd(sys.executable,spectra/"read_cubes.py",cache,"--positions",spectra/group,"--workers","2")
        for label,prior in receipts.items():
            got=json.loads((receipt_dir/(label+".json")).read_text())
            for key in ("tile","url","etag","size","sha256","start","end","pipeline"):
                assert got[key]==prior[key],f"{label}: changed {key}"
        cmd(sys.executable,spectra/"photometry.py",cache)
        SUMMARY["fresh_first_sample_tables"]=compare(old)
        SUMMARY["fresh_target_tiles"]=len(receipts)
        import csv as _csv
        _sys.path.insert(0,str(spectra))
        from read_cubes import fetch_tile
        from photometry import measure
        # Independent fresh data from the two distant, unresolved warm fields.
        with (STUDY/"fgk_spectra/positions.csv").open() as f: allpos=list(_csv.DictReader(f))
        selected=[r for r in allpos if r["label"].split("_")[0] in ("FGK18","FGK31")]
        assert len(selected)==6
        bytile={}
        for r in selected:bytile.setdefault(int(r["tile"]),[]).append(r)
        warmold={p.stem:json.loads(p.read_text()) for p in (STUDY/"fgk_spectra/cube_receipts").glob("*.json")}
        out=Path(tmp)/"warm_receipts";out.mkdir()
        warmcache=Path(tmp)/"warm";warmcache.mkdir()
        for tile,entries in bytile.items():
            fetch_tile((tile,entries,warmcache,out))
            got=json.loads((out/f"{tile}.json").read_text())
            prior=warmold[str(tile)]
            assert got["sha256"]==prior["sha256"] and got["size"]==prior["size"]
        spectra_rows,features,coverage,colors=measure([warmcache/(r["label"]+".npz") for r in selected])
        original_rows={ (r["label"],r["channel"]):r
                     for r in rows(STUDY/"fgk_spectra/channel_spectra.csv")}
        for r in spectra_rows:
            p=original_rows[r["label"],str(r["channel"])]
            for key in ("aperture_mJy","outlier_included_mJy"):
                a,b=value(r[key]),value(p[key])
                if a is None or not math.isfinite(a):
                    assert b is None or not math.isfinite(b)
                else:
                    assert b is not None and math.isclose(a,b,rel_tol=1e-5,abs_tol=1e-4), (r["label"],r["channel"],key,a,b)
        SUMMARY.update({"fresh_warm_tiles":len(bytile),"fresh_warm_positions":len(selected),
                        "fresh_warm_channels":len(spectra_rows)})
        # New byte-range read of the frozen blind tile and complete 58-star cohort.
        import importlib.util
        _sys.path.insert(0,str(STUDY/"blind"))
        import dr1_acquire
        blind=STUDY/"blind"
        before=json.loads((blind/"7119.json").read_text())
        with np.load(blind/"dr1_spectra_and_blanks.npz") as z:
            original_arrays={k:z[k].copy() for k in z.files}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            dr1_acquire.main()
        after=json.loads((blind/"7119.json").read_text())
        for field in ("sha256","etag","size","tile","url","start","end","pipeline"):
            assert before[field]==after[field],("blind",field)
        with np.load(blind/"dr1_spectra_and_blanks.npz") as z:
            assert set(z.files)==set(original_arrays)
            for k in z.files:
                np.testing.assert_allclose(z[k],original_arrays[k],atol=3e-5,rtol=1e-6,equal_nan=True,err_msg=k)
        SUMMARY["fresh_blind_tile"]=7119
        SUMMARY["fresh_blind_spectra"]=len(original_arrays)-1
    SUMMARY["status"]="passed"
    (STUDY/"referee_fresh.json").write_text(json.dumps(SUMMARY,indent=2,sort_keys=True)+"\n")
    print(json.dumps(SUMMARY,indent=2),flush=True)

if __name__=="__main__":
    if len(sys.argv)!=2 or sys.argv[1] not in ("offline","fresh"):
        raise SystemExit("usage: python referee_review.py [offline|fresh]")
    {"offline":offline,"fresh":fresh}[sys.argv[1]]()
