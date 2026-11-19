#!/usr/bin/env python3
"""Independent fail-closed replay and one authoritative status snapshot."""
from pathlib import Path
import csv,gzip,hashlib,json,math
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
BATCH=OUT/"beam_batches"
def sha(x):return hashlib.sha256(x).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+"\n")
def main():
    inp=OUT/"footprint_hits.csv"
    if not inp.exists():raise RuntimeError("No frozen footprint_hits.csv")
    with inp.open() as f:rows=list(csv.DictReader(f))
    batch_files=sorted(BATCH.glob("batch_????_manifest.json"))
    seen=set();recs=[];issues=[]
    for p in batch_files:
        meta=json.loads(p.read_text());f=BATCH/p.name.replace("_manifest.json",".jsonl.gz")
        if not f.exists() or sha(f.read_bytes())!=meta["data_sha256"]:
            issues.append("corrupt batch "+p.name);continue
        if meta["input_sha256"]!=sha(inp.read_bytes()):
            issues.append("input changed "+p.name);continue
        data=[json.loads(s) for s in gzip.decompress(f.read_bytes()).decode().splitlines() if s]
        if len(data)!=meta["n_positions"]:issues.append("bad count "+p.name)
        for offset,r in enumerate(data):
            idx=meta["first_row"]+offset
            if idx in seen:issues.append("duplicated input offset "+str(idx))
            seen.add(idx);recs.append(r)
    spec=sum(bool(r.get("desi_dr1",{}).get("matching_count")) for r in recs)
    photo=sum(bool(r.get("legacy_dr10_photoz",{}).get("matching_count")) for r in recs)
    hats=sum(any(x.get("status")=="OK" for x in r.get("euclid_hats",{}).get("catalogs",{}).values()) for r in recs)
    pixel=sum(bool(r.get("euclid_mer_beam_pixels",{}).get("full_four_band_science_pixels")) for r in recs)
    byband={b:sum(r.get("euclid_mer_beam_pixels",{}).get("bands",{}).get(b,{}).get("status")=="FULL_BEAM_NATIVE_PIXELS" for r in recs)
            for b in ("VIS","Y","J","H")}
    injections={b:sum(any(q.get("status")=="LOCAL_INJECTION_90PCT" for q in
               r.get("euclid_mer_beam_pixels",{}).get("bands",{}).get(b,{}).get(
                 "injections",{}).get("injections",[])) for r in recs) for b in byband}
    model=sum(isinstance(r.get("conditional_stellar_mass"),dict) and r["conditional_stellar_mass"].get("status")=="MODEL_CONDITIONAL_NOT_VALIDATED" for r in recs)
    matches=[r for r in recs if not r.get("desi_dr1",{}).get("matching_count") and
              not r.get("legacy_dr10_photoz",{}).get("matching_count") and
              r.get("euclid_mer_beam_pixels",{}).get("full_four_band_science_pixels")]
    failures={k:sum(r.get(k,{}).get("status")=="UNEXPECTED_STAGE_ERROR" for r in recs)
              for k in ("desi_dr1","legacy_dr10_photoz","euclid_hats","euclid_mer_beam_pixels")}
    complete=len(seen)==len(rows) and not issues
    vetted=sum(r.get("validated_starless",False) for r in recs)
    if vetted:issues.append("UNSUPPORTED VALIDATED STARLESS LABEL PRESENT")
    summary={"all_footprint_entries":len(rows),"analyzed_entries":len(recs),
       "unique_input_offsets":len(seen),"full_position_followup_complete":complete,
       "publication_ready":False,"disposition":"DRAFT_STILL_REQUIRED",
       "survey_counts":{"ALFALFA_alpha100":117,"FASHI_DR2":1258,"WALLABY_PDR2":0},
       "DESI_DR1_beam_velocity_matches":spec,"Legacy_DR10_photoz_credible_matches":photo,
       "Euclid_HATS_catalog_access_positions":hats,"full_four_band_beam_pixels":pixel,
       "full_beam_pixel_bands":byband,"local_injection_90pct_bands":injections,
       "model_conditional_Mstar_estimates_not_validated":model,
       "image_complete_no_matching_redshift_candidate_positions":len(matches),
       "image_complete_no_redshift_candidate_preliminary_list":[
           {"survey":x["survey"],"id":x["id"],"ra":x["ra_deg"],"dec":x["dec_deg"]}
           for x in matches[:100]],
       "validated_Mstar_upper_limits":0,"validated_MHI_Mstar_lower_bounds":0,
       "validated_starless_HI_galaxies":0,
       "fail_closed_reasons":[
           "MER science mosaic pixels lack independently verified per-exposure masks, PSF and zero-point propagation",
           "Blank-sky fixed-position exponential injections are not end-to-end all-beam blind-recovery completeness",
           "Single-dish HI centroid/linewidth may be tidally displaced or confused; cube/RFI/HVC vetting not supplied",
           "At HI redshift, photo-z membership is probabilistic and single-dish beam spectral association is not host identity",
           "No model-independent stellar population/dust/distance mass limit; conditional figures are not validated bounds",
           "Cloud-9/FASHI candidate aliases not completely reconciled, so no novelty claim"],
       "unexpected_stage_errors":failures,"integrity_warnings":issues}
    dump(OUT/"full_beam_referee.json",summary)
    note=("## All-position S3/HATS follow-up (independent status)\n\n"
       f"Processed **{len(recs)}/{len(rows)}** H I catalogue-entry beams in checkpoint batches "
       f"(complete: **{str(complete).lower()}**). DESI DR1 beam/velocity associations: **{spec}**; "
       f"Legacy DR10 photo-z credible-range associations: **{photo}**; "
       f"Euclid HATS accessed at **{hats}** positions; **{pixel}** have four-band native MER pixels "
       f"spanning one complete single-dish beam. Per-band full beam pixel counts: {byband}. "
       f"Local masked-injection 90% recovery estimates: {injections}. "
       f"Conditional mass scenarios: **{model}**, **none** validated. "
       f"Four-band pixel positions lacking these redshift associations: **{len(matches)}**, "
       "unvetted rather than proven starless.\n\n"
       "**Independent verdict:** zero validated starless galaxies, zero defensible Mstar upper limits "
       "or gas/stellar bounds; HI RFI, cubes, tidal environment, confusion, per-image quality flags "
       "and stellar SED priors prevent a physical discovery or a statistical incidence upper bound. "
       "**Remain draft; no [REMARKABLE] tag.** See "
       "[machine referee](output/full_beam_referee.json) and "
       "[per-position compressed batches](output/beam_batches/).\n")
    for p in (ROOT/"README.md",ROOT/"REFEREE.md"):
        content=p.read_text()
        begin="<!-- BEGIN ALL-BEAM FOLLOWUP -->"
        end="<!-- END ALL-BEAM FOLLOWUP -->"
        section=begin+"\n"+note+end
        if begin in content and end in content:
            a=content.index(begin);b=content.index(end,a)+len(end)
            content=content[:a]+section+content[b:]
        else:
            content=content.rstrip()+"\n\n"+section+"\n"
        p.write_text(content)
    root=ROOT.parent.parent/"README.md"
    if root.exists():
        s=root.read_text()
        token="- **Where Euclid's new images overlap clouds of hydrogen gas.**"
        for line in s.splitlines():
            if line.startswith(token):
                new=(token+f" The frozen HI survey inventory contains 1,375 beam/tile overlaps "
                   f"(117 ALFALFA, 1,258 FASHI, 0 WALLABY); {len(recs)}/1,375 received "
                   f"checkpointed S3/HATS/Legacy/DESI follow-up, with {pixel} full four-band-beam "
                   "MER pixel sets. A beam redshift association does not prove HI host identity; "
                   "no validated starless galaxy, stellar-mass limit or occurrence rate. "
                   "*Measured (incomplete physical vetting)* · "
                   "[Results and referee](results/euclid-hi-dark-census/README.md)")
                root.write_text(s.replace(line,new,1))
                break
    print("REFEREE",json.dumps({k:summary[k] for k in ("analyzed_entries","full_position_followup_complete",
        "DESI_DR1_beam_velocity_matches","full_four_band_beam_pixels","validated_starless_HI_galaxies")}),flush=True)
if __name__=="__main__":main()
