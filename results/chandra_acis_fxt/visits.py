"""Forced aperture photometry in independent released visits of retained leads.

python results/chandra_acis_fxt/visits.py --evidence results/chandra_acis_fxt/evidence
The prior blind selection is fixed. No trigger search, ancillary downloads, or
flux conversion is performed during this follow-up characterization.
"""
import argparse
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.stats import binom, chi2

from archive import dump, replay
from campaign import read, read_obs
from search import PIXEL_ARCSEC, _gti_map, _orientation, _read_events, _sky_wcs
from vet import separation


POINTING_RADIUS_ARCSEC = 5400.0
BIN_SECONDS = 60.0


def _lead_id(row):
    return f"{row['obsid']}/{row['input_sha256']}/{row['trigger_key']}"


def _receipt_plan(evidence, leads):
    """Filter on frozen detector reference positions, independent of source counts."""
    inventory = read(evidence/"inventory/inventory.json")
    jobs = []
    exclusions = []
    for obs in inventory:
        oid = int(obs["obsid"])
        path = evidence/"observations"/f"{oid}.json.gz"
        if not path.is_file():
            exclusions.append({"obsid":oid,"reason":"no observation ledger"})
            continue
        ledger = read_obs(path)
        if not ledger.get("acquired"):
            exclusions.append({"obsid":oid,"reason":"no verified EVT2 receipt"})
            continue
        for receipt in ledger["acquired"]:
            ra, dec = receipt.get("ra_nom_deg"), receipt.get("dec_nom_deg")
            if ra is None or dec is None or not np.isfinite([ra,dec]).all():
                exclusions.append({"obsid":oid,"sha256":receipt["sha256"],"reason":"missing FITS nominal pointing"})
                continue
            selected = [r for r in leads if int(r["obsid"]) != oid and
                        separation(ra,dec,r["ra_deg"],r["dec_deg"]) <= POINTING_RADIUS_ARCSEC]
            if selected:
                jobs.append((oid,receipt,selected))
    return inventory,jobs,exclusions


def _inside_segment(data, x, y, radius, a, b, duty, matrix, detector_limits):
    """Return counts over conservatively supported half-open bins of one CCD GTI."""
    ret = {"gti_met_s":[a,b], "duty":duty, "accepted_bins":[], "rejected_bins":{},
           "gross_counts":0, "annulus_counts":0}
    if matrix is None or len(data)<16 or b-a<30:
        ret["reason"]="insufficient events, good time or chip-to-sky geometry"
        return ret
    sx,sy=data["x"]-x,data["y"]-y
    dr=np.hypot(sx,sy)
    # Avoid selecting on a central aperture count. Independent photons must
    # sample the neighborhood around the aperture on all four sides.
    witness=(dr>4*radius)&(dr<max(150.,8*radius))
    quadrants=(sx[witness]>=0).astype(int)+2*(sy[witness]>=0).astype(int)
    quadrants=np.bincount(quadrants,minlength=4)
    ret["neighborhood_quadrant_counts"]=quadrants.tolist()
    if min(quadrants)<2:
        ret["reason"]="insufficient local footprint witnesses"
        return ret

    t=data["t"]
    edges=np.linspace(a,b,max(1,math.ceil((b-a)/BIN_SECONDS))+1)
    source = dr<=radius
    annulus = (dr>=2*radius)&(dr<=4*radius)
    on=np.histogram(t[source],edges)[0]
    off=np.histogram(t[annulus],edges)[0]
    far=np.flatnonzero(dr>4*radius)
    far=far[np.argsort(t[far],kind="stable")]
    times=t[far]
    # Each independent field event estimates where the fixed sky aperture
    # falls at its actual frame time. Use its full range, including dither.
    inverse=np.linalg.inv(matrix)
    center = (np.column_stack([data["cx"][far],data["cy"][far]]) -
              np.column_stack([sx[far],sy[far]]) @ inverse)
    cuts=np.searchsorted(times,edges)
    limits_low,limits_high=detector_limits
    guard=4*radius+4
    rejected=Counter()
    for i,(left,right) in enumerate(zip(edges,edges[1:])):
        begin,end=cuts[i:i+2]
        if end-begin<4 or np.diff(np.r_[left,times[begin:end],right]).max()>30:
            rejected["unsampled_dither"]+=1
            continue
        lo=center[begin:end].min(axis=0)-guard
        hi=center[begin:end].max(axis=0)+guard
        if (lo<limits_low).any() or (hi>limits_high).any():
            rejected["chip_or_subarray_edge"]+=1
            continue
        if any(lo[0]<=boundary<=hi[0] for boundary in (256.5,512.5,768.5)):
            rejected["node_boundary"]+=1
            continue
        ret["gross_counts"]+=int(on[i])
        ret["annulus_counts"]+=int(off[i])
        spans=ret["accepted_bins"]
        if spans and spans[-1][1]==float(left):
            spans[-1][1]=float(right)
        else:
            spans.append([float(left),float(right)])
    ret["rejected_bins"]=dict(rejected)
    if not ret["accepted_bins"]:
        ret["reason"]="no bins with supported on-chip aperture and annulus"
    return ret


def _measure(path, oid, leads):
    with fits.open(path,memmap=False) as hdul:
        evt=hdul["EVENTS"]
        h=evt.header
        if (int(h["OBS_ID"])!=oid or str(h["INSTRUME"]).strip().upper()!="ACIS" or
            str(h["GRATING"]).strip().upper()!="NONE" or str(h["READMODE"]).strip().upper()!="TIMED"):
            raise ValueError("Unexpected EVT2 instrument, readout, grating or ObsID")
        data=_read_events(evt)
        gtis=_gti_map(hdul)
        wcs=_sky_wcs(evt)
    valid_spans=[pair for spans in gtis.values() for pair in spans]
    if not valid_spans:
        raise ValueError("No time support in independent-visit event list")
    time_range=np.array([min(a for a,_ in valid_spans),max(b for _,b in valid_spans)])
    ref=h.get("MJDREF")
    if ref is None:
        ref=float(h["MJDREFI"])+float(h["MJDREFF"])
    epoch=Time(float(ref)+(time_range+float(h.get("TIMEZERO",0)))/86400,
               format="mjd",scale=h["TIMESYS"].strip().lower()).utc.mjd
    xy=wcs.all_world2pix([[r["ra_deg"],r["dec_deg"]] for r in leads],1)
    measurements={}
    for row,(x,y) in zip(leads,xy):
        # Keep one fixed sky aperture per source across all visits.
        radius=max(5.,2.*row["position_error_arcsec"])/PIXEL_ARCSEC
        relative=("before" if epoch[1]<=row["event_mjd"] else
                  "after" if epoch[0]>row["event_mjd"] else "overlapping_epoch")
        measurements[_lead_id(row)]={"sky_x":float(x),"sky_y":float(y),"radius_px":radius,
                                     "radius_arcsec":radius*PIXEL_ARCSEC,"segments":[],
                                     "event_mjd_utc":row["event_mjd"],"visit_mjd_utc":epoch.tolist(),
                                     "temporal_relation":relative}
    for ccd,intervals in gtis.items():
        chip=data[data["ccd"]==ccd]
        if not len(chip):
            continue
        low=np.array([chip["cx"].min(),chip["cy"].min()])
        high=np.array([chip["cx"].max(),chip["cy"].max()])
        matrix=_orientation(chip)
        on,live=h.get(f"ONTIME{ccd}"),h.get(f"LIVTIME{ccd}")
        duty=float(live)/float(on) if on and live is not None else float(h["DTCOR"])
        if not 0<duty<=1.00001:
            raise ValueError("Invalid CCD live fraction")
        for visit in measurements.values():
            x,y,r=visit["sky_x"],visit["sky_y"],visit["radius_px"]
            if not np.isfinite([x,y]).all():
                raise ValueError("Nonfinite independent-visit WCS projection")
            # A coordinate outside the empirical CCD bounding box cannot have
            # a supported source and annulus on all sides of the aperture.
            if not (chip["x"].min()<=x<=chip["x"].max() and chip["y"].min()<=y<=chip["y"].max()):
                continue
            for a,b in intervals:
                selected=chip[(chip["t"]>=a)&(chip["t"]<b)]
                piece=_inside_segment(selected,x,y,r,float(a),float(b),min(duty,1.),matrix,(low,high))
                piece["ccd_id"]=ccd
                visit["segments"].append(piece)
    for result in measurements.values():
        segments=result["segments"]
        spans=sorted((a,b,p["duty"]) for p in segments for a,b in p["accepted_bins"])
        if not spans:
            result["status"]="coverage_unconfirmed"
            continue
        if any(a<last_stop for (_,last_stop,_),(a,_,_) in zip(spans,spans[1:])):
            result["status"]="ambiguous_multiple_chip_coverage"
            continue
        live=sum((b-a)*d for a,b,d in spans)
        n=sum(p["gross_counts"] for p in segments)
        bg=sum(p["annulus_counts"] for p in segments)
        # Treat every aperture count as source for an upper limit conservative
        # with respect to additive background. The test of excess is separate.
        upper=float(.5*chi2.ppf(.95,2*(n+1)))
        result.update(status="measured_aperture",accepted_live_s=live,gross_counts=n,annulus_counts=bg,
                      annulus_area_ratio=12.,net_aperture_rate_ct_s=(n-bg/12)/live,
                      gross_counts_upper95=upper,conditional_aperture_rate_upper95_ct_s=upper/live,
                      local_excess_p=float(binom.sf(n-1,n+bg,1/13)))
    return measurements


def run_visits(evidence, leads, out):
    """Process each nearby byte-verified product once and reuse successful checkpoints."""
    if (os.getenv("GITHUB_ACTIONS")!="true" or
        os.getenv("GITHUB_REF")!="refs/heads/research/chandra-acis-fast-transients" or
        os.getenv("GITHUB_EVENT_NAME") not in {"push","workflow_dispatch"}):
        raise RuntimeError("Independent-visit retrieval requires the campaign workflow")
    evidence,out=Path(evidence),Path(out)
    out.mkdir(parents=True,exist_ok=True)
    if len({_lead_id(r) for r in leads})!=len(leads):
        raise ValueError("Lead identifiers must be unique")
    inventory,jobs,exceptions=(_receipt_plan(evidence,leads) if leads else ([],[],[]))
    deps=("visits.py","search.py","archive.py","vet.py","requirements.txt")
    hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in deps}
    code_hash=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()
    receipts=[]
    by_lead=defaultdict(list)
    for oid,receipt,candidates in jobs:
        identity=[code_hash,receipt["sha256"],[[_lead_id(r),r["ra_deg"],r["dec_deg"],
                                                r["position_error_arcsec"],r["event_mjd"]] for r in candidates]]
        key=hashlib.sha256(json.dumps(identity,separators=(",",":")).encode()).hexdigest()
        name=f"{oid}-{key[:20]}.json"
        checkpoint=out/name
        record=read(checkpoint) if checkpoint.is_file() else None
        if not record or record.get("status")!="complete" or record.get("key")!=key:
            record={"key":key,"obsid":oid,"input":receipt,"lead_ids":[_lead_id(r) for r in candidates]}
            try:
                with TemporaryDirectory() as directory:
                    tmp=Path(directory)
                    manifest=tmp/"one.json"
                    dump(manifest,[receipt])
                    replay(manifest,tmp)
                    record["measurements"]=_measure(tmp/receipt["filename"],oid,candidates)
                record["status"]="complete"
            except Exception as exc:
                record.update(status="error",error=f"{type(exc).__name__}: {exc}")
            dump(checkpoint,record)
        receipts.append({"obsid":oid,"input_sha256":receipt["sha256"],"record":name,"status":record["status"]})
        for candidate in candidates:
            index=_lead_id(candidate)
            m=record.get("measurements",{}).get(index,{})
            by_lead[index].append({"visit_obsid":oid,"record":name,"status":m.get("status",record["status"]),
                                  "temporal_relation":m.get("temporal_relation"),"local_excess_p":m.get("local_excess_p"),
                                  "conditional_aperture_rate_upper95_ct_s":
                                  m.get("conditional_aperture_rate_upper95_ct_s")})
    for visits in by_lead.values():
        for visit in visits:
            if visit["local_excess_p"] is not None:
                visit["per_lead_bonferroni_p"]=min(1.,visit["local_excess_p"]*len(visits))
    summary={"implementation_sha256":code_hash,"implementation_files_sha256":hashes,"leads_submitted":len(leads),"inventory_obsids":len(inventory),
             "pointing_prefilter_arcsec":POINTING_RADIUS_ARCSEC,"unassessable_receipts":exceptions,
             "nearby_products":receipts,"visits":{_lead_id(r):by_lead[_lead_id(r)] for r in leads},
             "interpretation":"Only other ObsIDs in the frozen release inventory enter this check. "
               "Positive time and count limits are conditional on a dither-sampled, geometrically "
               "interior aperture and annulus with independent neighboring field photons. Small bad "
               "pixels, exposure holes, vignetting and PSF aperture losses are uncalibrated. "
               "A count excess in a repeated pointing can also reflect diffuse or crowded emission. "
               "These limits apply to mean detected aperture counts during the stated covered bins."}
    dump(out/"independent_visits.json",summary)
    return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence",type=Path,required=True)
    args=parser.parse_args()
    run_visits(args.evidence,read(args.evidence/"leads.json"),args.evidence/"independent_visits")
