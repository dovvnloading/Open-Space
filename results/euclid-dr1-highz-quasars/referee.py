#!/usr/bin/env python3
"""Independent compact-data referee replay of frozen DR1 HATS cone selection.

Reads deposited broad rows, not import of the upstream selector; no web calls.
Verifies row accounting, flux cuts, signed VIS, control recovery, raw source
metadata SHA-256, publisher partition row counts, and output hash manifest.
Does NOT claim empirical purity, image fidelity, or spectroscopy.
"""
import csv,hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
D=ROOT/"s3_output"
def sha(data):return hashlib.sha256(data).hexdigest()
def val(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
def ab(v):return 23.9-2.5*math.log10(v) if v is not None and v>0 else None
def sep(a,b,c,d):
    x=math.radians(a-c);y=math.radians(b-d)
    s=math.sin(y/2)**2+math.cos(math.radians(b))*math.cos(math.radians(d))*math.sin(x/2)**2
    return 7200*math.degrees(math.asin(min(1,math.sqrt(max(0,s)))))
def try_independent(r):
    f={b:val(r.get("flux_"+b+"_3fwhm_aper")) for b in ("vis","y","j","h")}
    e={b:val(r.get("fluxerr_"+b+"_3fwhm_aper")) for b in f}
    sn={b:f[b]/e[b] if f[b] is not None and e[b] is not None and e[b]>0 else None for b in f}
    # NULL is never a nondetection. Signed VIS flux must have a valid error.
    if sn["vis"] is None or not -2<sn["vis"]<2:return False,None
    if sn["j"] is None or sn["j"]<6 or sn["h"] is None or sn["h"]<4:return False,None
    if r["flag_j"] not in ("0","0.0") or r["flag_h"] not in ("0","0.0"):return False,None
    el=val(r["ellipticity"])
    if el is None or el>=.3:return False,None
    ym=ab(f["y"]);jm=ab(f["j"]);hm=ab(f["h"])
    if jm is None or hm is None:return False,None
    jh=jm-hm
    if jh>.45:return False,jh
    if ym is None:
        if e["y"] is None or e["y"]<=0:return False,jh
        yj=ab(max(2*e["y"],f["y"] or -1))-jm
    else:yj=ym-jm
    if yj is None or yj<.40:return False,jh
    return True,jh

def audit():
    summary=json.loads((D/"summary.json").read_text())
    shortlist=json.loads((D/"shortlist.json").read_text())
    parts=json.loads((D/"partitions.json").read_text())
    freeze=json.loads((ROOT/"freeze.json").read_text())
    assert sha((ROOT/"freeze.json").read_bytes())==summary["original_freeze_sha256"]
    assert sha((ROOT/"run_s3.py").read_bytes())==summary["s3_script_sha256"]
    manifest=json.loads((D/"hashes.json").read_text())
    for p,m in manifest.items():
        path=ROOT/p
        assert path.is_file(),p
        data=path.read_bytes()
        assert len(data)==m["bytes"] and sha(data)==m["sha256"],p
    for kind in ("wide_survey","deep_survey"):
        content=(D/("partition_info_"+kind+".csv")).read_bytes()
        source=next(x for x in summary["source_metadata"] if x["kind"]==kind)
        assert sha(content)==source["sha256"]
        rows={(int(x["Norder"]),int(x["Npix"])):int(x["num_rows"]) for x in
              csv.DictReader(content.decode().splitlines())}
        for part in (p for p in parts if p["kind"]==kind and p["status"]=="READ"):
            key=(part["Norder"],part["Npix"])
            assert part["advertised_rows"]==rows[key]==part["num_rows"]
            assert part["s3_etag"] and part["object_uri"].endswith(".parquet")
    with (D/"broad_rows.csv").open(newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==summary["sources_satisfying_broad_NIR_VIStable_colours"]
    selection=[];dist=[]
    for row in rows:
        ok,jh=try_independent(row)
        if ok:
            selection.append(row)
            dist.append(jh)
    selection.sort(key=lambda x:(-float(x["flux_j_3fwhm_aper"])/
                                  float(x["fluxerr_j_3fwhm_aper"]),x["object_id"]))
    assert len(selection)==summary["strict_colour_shape_proxy_count"],(
        len(selection),summary["strict_colour_shape_proxy_count"])
    ids=[s["object_id"] for s in selection[:len(shortlist)]]
    assert ids==[s["id"] for s in shortlist],("Frozen shortlist priority mismatch",ids)
    control=next(x for x in summary["known_quasars"] if x["name"].startswith("EUCL J125308"))
    known=next(x for x in freeze["known"] if x["name"]==control["name"])
    near=[x for x in selection if sep(float(x["ra"]),float(x["dec"]),
                                      known["ra"],known["dec"])<1]
    assert bool(near)==control["dropout_cut_recovered"]
    assert summary["confirmed_new_quasars"]==0 and not summary["discovery_level"]
    requested=sum(2*math.pi*(180/math.pi)**2*(1-math.cos(math.radians(x["r"])))
                  for x in freeze["fields"])
    assert round(requested,6)==summary["nominal_requested_cones_deg2_NOT_VERIFIED_IMAGE_FOOTPRINT"]
    blue=sum(jh<-.5 for jh in dist)
    near_qso=sum(-.5<=jh<=.45 for jh in dist)
    report={
       "verdict":"REPRODUCIBLE_BOUNDED_COLOUR_SCREEN_NOT_QUASAR_DISCOVERY",
       "n_read_partitions":sum(p["status"]=="READ" for p in parts),
       "n_broad_replayed":len(rows),"n_strict_independently_replayed":len(selection),
       "n_extremely_blue_J_minus_H_lt_minus_0p5":blue,
       "n_plausible_JH_window_minus_0p5_to_0p45_POST_HOC_DIAGNOSTIC":near_qso,
       "n_known_z7p69_control_in_strict_pool":len(near),
       "control_angular_separation_arcsec":min(
           (sep(float(x["ra"]),float(x["dec"]),known["ra"],known["dec"]) for x in near),
           default=None),
       "shortlist_top5_ids_replayed":ids,
       "checks":["input freeze and script hashes","manifest hashes","HATS publisher metadata rows",
                 "independent signed VIS/NIR/flags/ellipticity recomputation",
                 "shortlist ordering","known quasar spatial recovery",
                 "unsearched area not an effective sky footprint"],
       "fatal_science_limits":["NIR roundness alone is not calibrated PSF morphology",
            "toy quasar/dwarf P has arbitrary prior and is NOT a posterior probability",
            "no real dwarf template library or red galaxy null calibration",
            "no multi-dither VIS artefact inspection",
            "unWISE W1/W2, proper motion, spectra, X-ray/radio and literature crossmatch NOT done",
            "five highest J-SNR candidates may be pathologically blue in J-H; no discovery claim"],
       "new_quasar_validated":False,"remarkable_label_permitted":False}
    (D/"referee.json").write_text(json.dumps(report,indent=2)+"\n")
    hashes={str(p.relative_to(ROOT)):{"bytes":p.stat().st_size,"sha256":sha(p.read_bytes())}
            for p in sorted(D.iterdir()) if p.is_file() and p.name!="hashes.json"}
    (D/"hashes.json").write_text(json.dumps(hashes,indent=2)+"\n")
    print("INDEPENDENT REFEREE",json.dumps({
      "broad":len(rows),"strict":len(selection),"extremely_blue":blue,
      "plausible_colour_window":near_qso,
      "known_control":bool(near),"verdict":report["verdict"]}))
if __name__=="__main__":audit()
