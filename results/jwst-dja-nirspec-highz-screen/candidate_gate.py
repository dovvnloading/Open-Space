#!/usr/bin/env python3
"""Fail-closed audit of a *human-reviewed* candidate evidence dossier.

This checks evidence completeness, NOT scientific correctness or source identity.
It does not generate positive dossiers or derive any 2D or NIRCam measurements.
Usage:
  python candidate_gate.py dossier.json
  python candidate_gate.py --self-test
Dossiers may only be prepared AFTER actual pixel-level processing; mere truthy
assertions without SHA-256 receipts will fail. Review the native pixels separately.
"""
import argparse
import json
import re
from pathlib import Path

SHA=re.compile(r"^[0-9a-f]{64}$")


def digest_ok(value):
    return isinstance(value,str) and bool(SHA.fullmatch(value))


def inspect(dossier):
    missing=[]
    source=dossier.get("source",{})
    for key in ("dja_1d_sha256","frozen_manifest_sha256"):
        if not digest_ok(source.get(key)):missing.append("source."+key)
    if not isinstance(source.get("spectrum_index"),int):
        missing.append("source.spectrum_index")
    if not (source.get("object_id") or source.get("coordinates_deg")):
        missing.append("source.object_id_or_coordinates")

    model=dossier.get("fit",{})
    for key in ("posterior_file_sha256","noise_calibration_sha256",
                "interloper_injection_receipt_sha256","stellar_templates_sha256"):
        if not digest_ok(model.get(key)):missing.append("fit."+key)
    try:
        lo=float(model["p_lowz_galaxy"])
        hi=float(model["p_z_gt_12"])
        z=float(model["z_median"])
        q16=float(model["z_q16"])
        q84=float(model["z_q84"])
        if not (0<=lo<=1 and 0<=hi<=1 and 0<=q16<=z<=q84):
            missing.append("fit.posterior_bounds")
        if model.get("posterior_calibrated") is not True:
            missing.append("fit.posterior_calibrated")
    except (KeyError,ValueError,TypeError):
        missing.append("fit.posterior_values")

    two=dossier.get("two_d",{})
    if not digest_ok(two.get("native_2d_sha256")):
        missing.append("two_d.native_2d_sha256")
    if two.get("lyman_break_visible_in_native_2d") is not True:
        missing.append("two_d.lyman_break_visible_in_native_2d")
    visits=two.get("independent_exposures",[])
    if len(visits)<2 or len({v.get("exposure_id") for v in visits if isinstance(v,dict)})<2:
        missing.append("two_d.independent_exposures")
    for index,v in enumerate(visits):
        if not isinstance(v,dict) or not digest_ok(v.get("file_sha256")):
            missing.append(f"two_d.independent_exposures[{index}].file_sha256")
        elif v.get("break_consistent") is not True:
            missing.append(f"two_d.independent_exposures[{index}].break_consistent")
    if two.get("bad_pixels_reviewed") is not True:
        missing.append("two_d.bad_pixels_reviewed")
    if two.get("overlapping_shutters_reviewed") is not True:
        missing.append("two_d.overlapping_shutters_reviewed")
    if two.get("unresolved_contamination") is not False:
        missing.append("two_d.unresolved_contamination")

    phot=dossier.get("nircam",{})
    if not digest_ok(phot.get("catalog_flux_limits_sha256")):
        missing.append("nircam.catalog_flux_limits_sha256")
    if not digest_ok(phot.get("filter_curve_sha256")):
        missing.append("nircam.filter_curve_sha256")
    if phot.get("highz_sed_consistent") is not True:
        missing.append("nircam.highz_sed_consistent")
    if phot.get("lowz_alternatives_fitted") is not True:
        missing.append("nircam.lowz_alternatives_fitted")

    lit=dossier.get("priority",{})
    if not digest_ok(lit.get("search_receipt_sha256")):
        missing.append("priority.search_receipt_sha256")
    if not all(lit.get(key) is True for key in
               ("ads","arxiv","dja","jades","ceers","uncover","mom","glass")):
        missing.append("priority.all_required_sources_searched")
    if lit.get("identity_classification") not in ("known","new_candidate","uncertain"):
        missing.append("priority.identity_classification")

    # The gate only permits further independent review; it never certifies a
    # discovery, even with a fully populated dossier.
    return {"status":("DOCUMENTATION_PRESENT_REQUIRES_INDEPENDENT_REFEREE"
                      if not missing else "INCOMPLETE_NOT_DISCOVERY"),
            "missing":missing,"passes_required_document_fields":not missing,
            "not_established":["redshift physical correctness","absence of slit contamination",
                               "photometric identity","novelty","discovery status"]}


def self_test():
    none=inspect({})
    assert none["status"]=="INCOMPLETE_NOT_DISCOVERY"
    assert len(none["missing"])>=10
    false=inspect({"two_d":{"lyman_break_visible_in_native_2d":True},
                   "fit":{"p_lowz_galaxy":0,"p_z_gt_12":1,
                          "z_median":15,"z_q16":14.5,"z_q84":15.5}})
    assert false["status"]=="INCOMPLETE_NOT_DISCOVERY"
    print("PASS: empty/unsupported high-z dossiers cannot pass the gate")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("dossier",nargs="?")
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:self_test()
    else:
        if not args.dossier:p.error("provide a candidate JSON dossier")
        result=inspect(json.loads(Path(args.dossier).read_text()))
        print(json.dumps(result,indent=2))
        if result["status"]=="INCOMPLETE_NOT_DISCOVERY":
            raise SystemExit(2)
