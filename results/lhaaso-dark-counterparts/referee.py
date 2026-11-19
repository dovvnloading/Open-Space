"""Independent in-repository referee for the 1LHAASO public-catalogue screening."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
CAT=OUT/"catalogs"

def main():
    receipts=json.loads((OUT/"receipts.json").read_text())
    summary=json.loads((OUT/"summary.json").read_text())
    checks=[]
    def check(name,ok,detail):
        checks.append({"name":name,"pass":bool(ok),"detail":detail})
    primary_receipt=receipts["lhaaso"]
    frozen=OUT/primary_receipt.get("frozen_path","frozen/1lhaaso_cao2024_vizier.csv")
    sha=hashlib.sha256(frozen.read_bytes()).hexdigest()
    check("original-bytes-sha256",sha==receipts["lhaaso"]["raw_sha256"],sha)
    is_preprint=primary_receipt.get("input_variant")=="gammapy_2023_preprint_transcription"
    primary=pd.read_csv(frozen,comment=None if is_preprint else "#",low_memory=False)
    unique_names=(primary.iloc[:,0].astype("string").str.strip().replace("",pd.NA).ffill().nunique() if is_preprint
                  else primary["1LHAASO"].nunique())
    source=pd.read_csv(CAT/"1lhaaso_sources.csv")
    check("180-components-90-sources",len(primary)==180 and
          unique_names==90 and len(source)==90,
          {"input_rows":len(primary),"input_sources":unique_names,
           "variant":primary_receipt.get("input_variant")})
    check("43-uhe",int(source.uhe.sum())==43,{"uhe":int(source.uhe.sum())})
    check("unique-position-valid",source.name.is_unique and source.ra.between(0,360).all() and
          source.dec.between(-90,90).all(),"unit: published source, not detector")
    dark=source.loc[~source.known_tev]
    check("dark-sample-recount",len(dark)==summary["published_tev_unassociated"],len(dark))
    class_rows=pd.read_csv(OUT/"source_class_controls.csv")
    n_classes=class_rows["class"].nunique()
    check("all-dark-sources-screened-or-marked-pending",n_classes==9 and len(class_rows)==len(dark)*n_classes,
          {"found":len(class_rows),"expected":len(dark)*n_classes,"classes":n_classes})
    active=class_rows[class_rows.status.str.startswith("screened")].copy()
    check("some-real-class-coverage",len(active)>0,{"eligible_rows":len(active)})
    if len(active):
        n=pd.to_numeric(active.control_n,errors="coerce")
        p=pd.to_numeric(active.chance_any_in_radius,errors="coerce")
        q=pd.to_numeric(active.p_nearest_margin,errors="coerce")
        check("finite-sample-p-bounds",(n.between(64,512)&p.ge(1/(n+1)-1e-9)&
              p.le(1)&q.ge(1/(n+1)-1e-9)&q.le(1)).all(),
              {"n_min":int(n.min()),"n_max":int(n.max())})
    erass=class_rows[class_rows["class"]=="eRASS1_hard_Xray"]
    if len(erass):
        out=erass[erass.status=="outside_survey_footprint"]
        check("unobserved-erass-not-nondetection",not out.n_matched.notna().any(),
              {"out_of_footprint":len(out)})
    strict=pd.read_csv(OUT/"strict_dark.csv")
    incomplete=strict[strict.status!="screen_complete"]
    check("strict-dark-is-tristate",not incomplete.no_geometrical_GeV_PSR_SNR.notna().any(),
          {"incomplete":len(incomplete)})
    pairpath=OUT/"candidate_pairs.csv"
    if pairpath.stat().st_size>10:
        pairs=pd.read_csv(pairpath)
        if len(pairs):
            check("pair-within-morphology-aperture",
                  (pairs.separation_deg<=pairs.source_screen_radius_deg+
                   pairs.counterpart_radius_deg+1e-5).all(),{"pairs":len(pairs)})
            check("no-pair-physical-claim",
                  pairs.association_status.eq("geometrical_only").all(),
                  "All candidate labels geometrical")
    missing=sorted(k for k,v in receipts.items() if v.get("status")!="loaded")
    check("no-novel-class-claim",not summary["discovery_claim"] and
          summary["physical_identifications"]==0,"Claim gate unchanged")
    result={"checks":checks,"pass":all(x["pass"] for x in checks),
            "passed":sum(x["pass"] for x in checks),"total":len(checks),
            "missing_catalogues":missing,"all_requested_catalogues_complete":not bool(missing),
            "discovery_level_finding":False,"remarkable_permitted":False,
            "independence":"independent code/data bookkeeping review, not independent telescope measurement"}
    (OUT/"referee.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    if not result["pass"]:
        raise SystemExit("Referee assertion failed: "+repr([x for x in checks if not x["pass"]]))
    census=pd.read_csv(OUT/"class_census.csv").fillna("")
    lit=pd.read_csv(ROOT/"literature.csv")
    rows=[
      "# 1LHAASO counterpart census: checked run output",
      "",
      "**Result:** positional/energetic candidate screen only; no new accelerator identified.",
      "",
      "Frozen Cao 1LHAASO transcription: "+str(len(source))+" physical sources, "+
      str(int(source.uhe.sum()))+" UHE, "+str(len(dark))+
      " lacking an original TeV association; "+str(int(dark.uhe.sum()))+" of the latter UHE.",
      "A negative association is not a survey non-detection unless that class has a validated footprint.",
      ("**Primary-data limitation:** Gammapy transcription of the 2023 Cao preprint is frozen and hashed; "
       "the 2024 CDS VizieR bytes could not be verified during the Nov 2026 CDS outage."
       if is_preprint else "The 2024 CDS VizieR original mirror bytes are archived and hashed."),
      "",
      "## Sky-shifted class controls",
      "",
      "| Class | Status | Covered sources | Observed matches | Mean control matches | Tail p |",
      "|---|---|---:|---:|---:|---:|",
    ]
    for _,r in census.iterrows():
        rows.append("| "+" | ".join(str(r.get(k,"")) for k in
             ["class","status","eligible_sources","observed_matches",
              "null_mean_matches","p_global_ge_observed"])+" |")
    rows += [
      "",
      "Individual source/class probabilities: source_class_controls.csv. "
      "Pair offsets, conditional chance rates and pulsar efficiency/size tests: candidate_pairs.csv.",
      "Chance probability uses (1 + matched sky-shift controls) / (1 + eligible controls). "
      "512 shared longitude offsets; Galactic latitude held fixed, shifts >= 5 degrees. "
      "Source aperture = clip(ePos95stat + r39 + 0.2 deg, 0.35 deg, 1.5 deg); "
      "catalogued counterpart radius added.",
      "",
      "## Critical limitations",
      "",
      "- This chance test matches Galactic latitude but not local longitude, source confusion "
      "or exposure in each survey. It is a screening null, not a source-identity posterior.",
      "- eROSITA/SMGPS/HESS survey masks are approximate; their global p-values are withheld. "
      "eRASS1 hard is not the 930203-source soft main catalogue. Only the German half-sky is public.",
      "- SMGPS SNR candidates are not the 510599-source compact radio catalogue. "
      "The project's prior SMGPS image-index astrometry correction is relevant when using maps.",
      "- 4HWC (Alfaro et al. 2026 arXiv:2602.00263) is the newest HAWC catalogue; "
      "its TeV detections and association proposals pre-empt novelty. "
      "Declination acceptance is approximate; no global HAWC4 p-value is quoted.",
      "- Fermi DR4 v35 is a corrected 2024 release of a catalogue first published 2023; "
      "HGPS is 2018 legacy. Never attribute new detections to a survey they predate.",
      "- ATNF spin-down and distance allow a necessary power/size test, not causality; "
      "KM2A 25-100 TeV power-law flux is not a measured bolometric budget.",
      "- Multiple correlated source/class hypotheses and crowded Galactic longitude "
      "invalidate a naive discovery significance from a small nominal p.",
      "",
      "## Pre-existing 2024-2026 literature claims (finite, not exhaustive)",
      "",
      "| 1LHAASO source | Reference | Earlier suggestion |",
      "|---|---|---|",
    ]
    for _,r in lit.iterrows():
        rows.append("| "+str(r["source"])+" | [paper]("+str(r["url"])+") | "+
                    str(r["previous_proposal"])+" |")
    rows += ["",
       "## Independent referee",
       "",
       str(result["passed"])+"/"+str(result["total"])+" assertions passed. "
       "Discovery-level label: **not permitted**.",
       "Pending catalogue/coverage checks: "+(", ".join(missing) if missing else "none registered")+".",
       "Still requires independent TeV morphology/spectrum, all-region exposure masks, "
       "full per-source literature checks, and targeted multiwavelength observations.",
       "",
       "Frozen source bytes, normalized catalogue extracts, byte hashes and original URLs "
       "are preserved in this folder; see receipts.json.",
       ""]
    (OUT/"REPORT.md").write_text("\n".join(rows))
    print("REFEREE PASS",result["passed"],"/",result["total"],"CLAIM GATE: NOT REMARKABLE")

if __name__=="__main__":
    main()
