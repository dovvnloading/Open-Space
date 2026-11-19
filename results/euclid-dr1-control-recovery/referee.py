#!/usr/bin/env python3
"""Independent conservative referee: original source cuts, frozen control
selection via scalar arithmetic, all reselected partition receipts and all
ranked candidates. No self-certification of spectroscopy or image clearance.
"""
import gzip,hashlib,json,math
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
OLD=ROOT.parent/"euclid-dr1-wide-highz-vetting"
def jd(p):return json.loads(p.read_text())
def hashb(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scalar(r):
    s=r["source"]
    def f(k):
        try:
            x=float(s.get(k))
            return x if math.isfinite(x) else None
        except (ValueError,TypeError):return None
    def sn(b):
        a,e=f("flux_"+b+"_3fwhm_aper"),f("fluxerr_"+b+"_3fwhm_aper")
        return a/e if a is not None and e is not None and e>0 else None
    j,h,y=f("flux_j_3fwhm_aper"),f("flux_h_3fwhm_aper"),f("flux_y_3fwhm_aper")
    ye=f("fluxerr_y_3fwhm_aper")
    v,jj,hh=sn("vis"),sn("j"),sn("h")
    if j is None or j<=0 or h is None or h<=0:return False
    col=(-2.5*math.log10(y/j) if y is not None and y>0 else
         -2.5*math.log10(2*ye/j) if ye is not None and ye>0 else -math.inf)
    jh=-2.5*math.log10(j/h)
    spur=f("spurious_prob"); el=f("ellipticity"); mu=f("mumax_minus_mag")
    return bool(j>.65 and hh is not None and hh>=4 and
        jj is not None and jj>=6 and v is not None and abs(v)<7 and
        col>=0 and -.20<=jh<=.48 and el is not None and el<.4 and
        mu is not None and -1.45<=mu<=-.75 and (spur is None or spur<.7))

def main():
    report={"status":"FAIL_UNTIL_EVIDENCE","discovery_level":False,"checks":[],"warnings":[]}
    q=jd(OUT/"control_audit.json")
    control=jd(OUT/"calibration.json")
    frozen=jd(ROOT/"freeze.json")
    matched=[r for r in q if "source" in r]
    matchedq=[r for r in matched if r["class"]=="quasar"]
    matchedd=[r for r in matched if r["class"]=="dwarf"]
    if (len(matchedq),len(matchedd))!=(24,14):raise AssertionError("Changed control sample")
    if sum(r["original_strict_pass"] for r in matchedq)!=3:
        raise AssertionError("Original controls did not replay")
    qpass=sum(scalar(r) for r in matchedq);dpass=sum(scalar(r) for r in matchedd)
    if (qpass,dpass)!=(21,4):
        raise AssertionError("Frozen scalar recovery not 21/24,5/14: "+str((qpass,dpass)))
    if qpass!=control["matched_quasar_recovery"]["numerator"]:
        raise AssertionError("Candidate selector and referee differ")
    report["checks"].append("Scalar independent control replay: 21/24 matched quasars, 4/14 matched dwarfs; old 3/24")
    seen=set();counts=Counter()
    source_etags={}
    oldsum=jd(OLD/"output/summary.json")
    for b in range(32):
        old=jd(OLD/f"output/batches/batch_{b:03d}_manifest.json")
        for p in old["partitions"]:
            source_etags[(p["Norder"],p["Npix"])]=(p["key"],p["etag"],p["rows"])
    if len(source_etags)!=768:raise AssertionError("Original source partition list damaged")
    for b in range(32):
        m=jd(OUT/f"batches/batch_{b:03d}_manifest.json")
        gz=OUT/f"batches/batch_{b:03d}_top.jsonl.gz"
        if m["top_sha256"]!=hashb(gz):raise AssertionError("Corrupt source batch "+str(b))
        if m["selection_id"]!=frozen["selection_id"]:raise AssertionError("Version drift")
        with gzip.open(gz,"rt",encoding="utf8") as f:
            rows=[json.loads(line) for line in f]
        if len(rows)!=m["retained_top_rows"]:raise AssertionError("Batch size mismatch")
        by=Counter((r["Norder"],r["Npix"]) for r in rows)
        if any(n>200 for n in by.values()):raise AssertionError("Partition top cap")
        for p in m["partitions"]:
            item=(p["Norder"],p["Npix"])
            if item in seen:raise AssertionError("Duplicate partition")
            seen.add(item)
            if (p["key"],p["etag"],p["read_rows"])!=source_etags[item]:
                raise AssertionError("Publisher identity/rows inconsistent")
            if by[item]!=p["retained_top_rows"]:
                raise AssertionError("Retained row count mismatch")
            counts["rows"]+=p["read_rows"];counts["eligible"]+=p["eligible_rows"]
            counts["ranked"]+=p["retained_top_rows"]
    if len(seen)!=768 or counts["rows"]!=888152165:
        raise AssertionError("Incomplete source rescan")
    summary=jd(OUT/"summary.json")
    if counts["eligible"]!=summary["eligible_catalogue_rows"] or summary["scanned_partitions"]!=768:
        raise AssertionError("Summary mismatched raw batch receipts")
    ranked=jd(OUT/"top_candidates.json")
    if len(ranked)!=100 or len({r["id"] for r in ranked})!=100:
        raise AssertionError("Truncated top 100")
    if any(r.get("spectroscopically_new_Q") or r.get("image_artefacts_cleared") for r in ranked):
        raise AssertionError("Unsupported discovery or image clearance")
    report["checks"].append("32 frozen, hash-verified HATS source batches: 768/768; 888152165/888152165 publisher rows")
    report["checks"].append("All selected counts, rank-preserving bounded top-200-per-partition registers reconcile")
    pixels=jd(OUT/"image_pixel_receipts.json") if (OUT/"image_pixel_receipts.json").exists() else []
    external=jd(OUT/"external_summary.json") if (OUT/"external_summary.json").exists() else {}
    goodpix=sum(any(v.get("status")=="REAL_FITS_SCIENCE_PIXELS_RANGE_READ_NOT_HUMAN_CLEARED"
             and v.get("science_pixels",0)>40 and v.get("science_pixel_sha256")
             for v in x.get("bands",[])) for x in pixels)
    if any(x.get("image_artefacts_cleared") for x in pixels):
        raise AssertionError("Pixel receipt falsely declares artefact clearance")
    report["actual_FITS_pixel_targets"]=goodpix
    report["image_queries"]=len(pixels)
    if goodpix<8:report["warnings"].append("Not all eight top ranks obtained usable science FITS pixels; artefacts remain unresolved")
    if external.get("actual_public_spectra_examined",0)==0:
        report["warnings"].append("Public spectrum *metadata* only; no new redshift established")
    model=jd(OUT/"conditional_model.json")
    two_class=jd(OUT/"ranked_conditional_top100.json")
    if len(two_class)!=100 or model["field_quasar_probability_calibrated"] is not False:
        raise AssertionError("Conditional model incorrectly reported field calibration")
    if model["training_label_quasars"]!=24 or model["training_label_dwarfs"]!=14:
        raise AssertionError("Unexpected label cohort")
    pcount=0
    for i,r in enumerate(two_class):
        p=r.get("conditional_p_quasar_vs_dwarf")
        pd=r.get("conditional_p_dwarf_vs_quasar")
        if p is not None:
            pcount+=1
            if pd is None or abs((p+pd)-1)>3e-8:
                raise AssertionError("Conditional two-class posterior fails to sum to one")
            if r["conditional_colour_density_status"]!="UNVALIDATED_CONDITIONAL_TWO_CLASS_MODEL":
                raise AssertionError("Unsupported extrapolated posterior")
        elif r["conditional_colour_density_status"]=="UNVALIDATED_CONDITIONAL_TWO_CLASS_MODEL":
            raise AssertionError("Claimed numeric model without posterior")
        if r.get("field_quasar_probability") is not None:
            raise AssertionError("Invented field prior")
    if pcount!=model["top100_with_numeric_conditional_probability"]:
        raise AssertionError("Conditional-model scored count mismatch")
    report["conditional_two_class_model_n_supported"]=pcount
    report["conditional_model_LOO_Brier"]=model["cross_validated_leave_one_out_Brier"]
    report["checks"].append("Top-100 quasar-versus-dwarf probabilities are explicitly conditional and withheld outside labelled colours")
    mirror=jd(OUT/"wise_mirror_summary.json") if (OUT/"wise_mirror_summary.json").exists() else None
    if mirror is not None:
        report["CatWISE_mirror_matches"]=mirror["CatWISE2020_mirror_position_matched_W1_W2"]
        report["unWISE_mirror_matches"]=mirror["unWISE_fallback_position_matches"]
    report["warnings"]+= [
      "Training and evaluation use same labelled control cohort: 21/24 is conditional, resubstitution apparent recall, not calibrated field completeness",
      "4/14 known dwarfs retained; weak and nonrepresentative contamination estimate; galaxies/blends/other dwarfs absent",
      "Catalogue HATS cells are not a clean imaging effective-area map",
      "No spectroscopic new quasar. A photometric ranked source is not a discovery."]
    report.update(status="INDEPENDENT_RECEIPT_PASS_PHOTOMETRIC_ONLY_NO_DISCOVERY",
      selection_id=frozen["selection_id"],
      checked_source_rows=counts["rows"],selected_rows=counts["eligible"],
      known_matched_quasars_recovered=qpass,matched_dwarfs_retained=dpass,
      field_completeness_calibrated=False,
      new_spectroscopic_quasar_confirmations=0,
      first_prior_PR_208_read_only=True)
    (OUT/"referee.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("REFEREE",json.dumps(report),flush=True)
if __name__=="__main__":main()
