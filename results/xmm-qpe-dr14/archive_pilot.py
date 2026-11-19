#!/usr/bin/env python3
"""Deterministic three-obs pilot of frozen QPE finder; NOT complete survey or gate pass."""
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile

from astropy.io import fits
import control_validation as cv

BASE=cv.BASE
OUT=cv.OUT
cv.MAX_SECONDS_PER_TARGET=5*60
cv.MAX_TAR_BYTES=3*1024**3
TARGETS_PER_JOB=3
POSITIONS=OUT/"archive_pn_point_freeze.csv"
FROZEN=OUT/"archive_pn_point_freeze.json"

def verified_manifest():
    frozen=json.loads(FROZEN.read_text())
    assert frozen["status"]=="PROVISIONAL_SAMPLE_FROZEN"
    assert cv.checksum(POSITIONS)==frozen["manifest_sha256"]
    rows=list(csv.DictReader(POSITIONS.open(newline="")))
    assert len(rows)==frozen["sdss_nuclear_point_pn_detections"]==1108
    return rows,frozen

def plan(rows):
    # Choose BEFORE accessing PPS: point-like, pn >=1000 photons,
    # offset<=1arcsec, catalogue pile-up diagnostic below 0.05.
    candidates=[r for r in rows if
                float(r["pn_pileup_diagnostic"] or 1)<.05 and
                float(r["pn_8_cts"])>=1000 and
                float(r["offset_arcsec"])<=1]
    candidates.sort(key=lambda r:(-float(r["score"]),
                                  -float(r["pn_8_cts"]),r["detid"]))
    seen=set();selected=[]
    for r in candidates:
        if r["obs_id"] not in seen:
            selected.append(r)
            seen.add(r["obs_id"])
        if len(selected)==TARGETS_PER_JOB:break
    assert len(selected)==TARGETS_PER_JOB
    return selected,{"rule":"Top catalogue variability rank, PN photons and clean low-pile-up diagnostic; no PPS data looked at",
        "parent_eligible_detections":len(rows),"ranking_preselection_count":len(candidates),
        "sample_obsids":[r["obs_id"] for r in selected],
        "primary_rows":[{k:r[k] for k in ("detid","obs_id","srcid","score",
                             "pn_8_cts","pn_pileup_diagnostic","offset_arcsec","ra","dec")}
                        for r in selected],
        "source_manifest_sha256":cv.checksum(POSITIONS),
        "scientific_status":"PILOT ONLY; cannot infer archival prevalence, FPR or null"}

def parse_obs(tarpath,obsid,rows):
    output={r["detid"]:{"detid":r["detid"],"obs_id":obsid,
         "ra":r["ra"],"dec":r["dec"],"srcid":r["srcid"],
         "status":"no_verified_PN_NET_curve",
         "pn_pileup_catalogue_diagnostic":r["pn_pileup_diagnostic"]}
        for r in rows}
    with tarfile.open(tarpath,"r:*") as tar:
        members=[m for m in tar if m.isfile() and
          Path(m.name).name.upper().startswith("P"+obsid+"PNX000SRCTSR8")
          and Path(m.name).suffix.upper() in (".FTZ",".FIT",".FITS")
          and m.size<=100*1024**2]
        for m in members:
            try:
                name=Path(m.name).name
                raw=tar.extractfile(m).read()
                if name.upper().endswith(".FTZ"):raw=gzip.decompress(raw)
                with fits.open(io.BytesIO(raw),memmap=False) as hdul:
                    data=next((x for x in hdul[1:] if x.data is not None
                         and hasattr(x.data,"names") and "RATE1" in x.data.names),None)
                    if data is None:continue
                    pos=cv.source_sky_position([data.header,hdul[0].header])
                    if not pos:continue
                    nearby=sorted((cv.sep_arcsec(pos,(float(r["ra"]),float(r["dec"]))),r)
                                  for r in rows)
                    sep,best=nearby[0]
                    if sep>15 or (len(nearby)>1 and nearby[1][0]-sep<2):
                        continue
                    if output[best["detid"]]["status"]!="no_verified_PN_NET_curve":
                        continue
                    if str(data.header.get("HDUCLAS2","")).strip()!="NET":
                        output[best["detid"]]["status"]="not_NET";continue
                    t,y,e,b,f,hard,td=cv.combine_energies(data)
                    curve=cv.rebin(t,y,e,b,f,hard,td)
                    match=cv.find_eruptions(curve,cv.mjd_reference(data.header),(.5,24.0))
                    row=output[best["detid"]]
                    row.update(status="source_position_verified_pn_net",
                        fits_source_position_deg=list(pos[:2]),
                        separation_arcsec=round(sep,3),
                        product=name,product_sha256=hashlib.sha256(raw).hexdigest(),
                        n_good_bins=int(curve["good"].sum()),
                        n_eruptions=match.get("n_peaks"),
                        peak_finder=match,
                        recurrence_unvetted=bool(match.get("range_match") and
                                                  match.get("n_peaks",0)>=2),
                        independent_particle_background_veto="NOT PERFORMED",
                        pileup_epatplot="NOT PERFORMED",
                        stellar_blazar_literature_vetting="NOT PERFORMED")
            except Exception as exc:
                # Another source product may be valid; never classify a parse failure as negative.
                output.setdefault("__errors__",{"errors":[]})["errors"].append({
                    "file":Path(m.name).name,"error":repr(exc)[:200]})
    return list(output.values()),len(members)

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg.get("stage")=="archive_pilot"
    rows,frozen=verified_manifest()
    primaries,selection=plan(rows)
    OUT.joinpath("archive_pilot_plan_01.json").write_text(
        json.dumps(selection,indent=2)+"\n")
    targets=selection["sample_obsids"]
    receipt={"as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "status":"running","pilot_obsids":targets,
        "limits":"5 min per obs, 3 obs per job, 3GiB per obs; downloaded tar and source product SHA256",
        "n_eligible_obsids":frozen["unique_obsids"],
        "parent_manifest_sha256":frozen["manifest_sha256"],
        "method":"FROZEN NET EPIC-pn 300s control_validation.find_eruptions",
        "processed_sources":[],"downloads":[]}
    for obsid in targets:
        local=[r for r in rows if r["obs_id"]==obsid]
        print("PILOT OBS",obsid,"eligible sources",len(local),flush=True)
        with tempfile.TemporaryDirectory() as td:
            tar=cv.download_product(obsid,"pilot:"+obsid,Path(td))
            if tar is None:
                receipt["processed_sources"].extend(
                    [{"detid":r["detid"],"obs_id":obsid,"status":"PPS_unavailable"} for r in local])
                receipt["downloads"].append(cv.RECEIPTS[-1])
            else:
                try:
                    data,nproducts=parse_obs(tar,obsid,local)
                    receipt["processed_sources"].extend(data)
                    receipt["downloads"].append({**cv.RECEIPTS[-1],
                                                   "n_merged_pn_products":nproducts})
                except Exception as ex:
                    receipt["processed_sources"].extend(
                        [{"detid":r["detid"],"obs_id":obsid,
                          "status":"PPS_parse_error","error":repr(ex)[:200]} for r in local])
                    receipt["downloads"].append(cv.RECEIPTS[-1])
        OUT.joinpath("archive_pilot_batch_01.json").write_text(
            json.dumps(receipt,indent=2,allow_nan=False,default=str)+"\n")
    valid=[r for r in receipt["processed_sources"]
           if r["status"]=="source_position_verified_pn_net"]
    flags=[r for r in valid if r.get("recurrence_unvetted")]
    receipt.update(status="PILOT_FINISHED_NOT_GATED",eligible_obsids=len(targets),
        obsids_with_parsed_pn=len({r["obs_id"] for r in valid}),
        pn_detections_screened=len(valid),
        unvetted_recurrent_flags=len(flags),validated_new_qpe=0,
        survey_incomplete_eligible_obsids=frozen["unique_obsids"]-len({r["obs_id"] for r in valid}),
        archive_screen_permitted=False,
        limitations=["No independent 4XMM-host LS all-sky census; Gaia/blazar exclusions not complete",
                     "No independent background region, particle GTI, epatplot pile-up veto, or spectral evolution",
                     "Only first three pre-registered catalogue-ranked observations attempted",
                     "No archive-wide empirical false-positive rate or injection-derived selection efficiency",
                     "Flagged source positions are not validated new QPEs; TDE/AGN/stellar novelty untested"])
    OUT.joinpath("archive_pilot_batch_01.json").write_text(
         json.dumps(receipt,indent=2,allow_nan=False,default=str)+"\n")
    print(json.dumps({k:receipt[k] for k in ("status","pn_detections_screened",
          "unvetted_recurrent_flags","obsids_with_parsed_pn",
          "survey_incomplete_eligible_obsids")},indent=2),flush=True)

if __name__=="__main__":
    main()
