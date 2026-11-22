"""Inspect ACIS readout frames for an existing, byte-pinned primary real trigger."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from zipfile import ZipFile

import numpy as np
from astropy.io import fits
import campaign
import diagnostics
import followup
import ingest_primary as ingest
import upload_recovery

OID = 24912
EVENT = "7be4b001269694c0680f601ae3750a54be83a9fdac18f01fbbdcecf32b6966fd"
TAG = "real/0"
GROUP = 91
need,sha,save = ingest.require,ingest.sha,ingest.save
DOC = "https://cxc.harvard.edu/ciao/data_products_guide/events.html"


def digest(b):
    return hashlib.sha256(b).hexdigest()


def csvz(path, names, rows):
    with path.open("wb") as out, gzip.GzipFile(filename="",mode="wb",fileobj=out,mtime=0) as packed:
        with io.TextIOWrapper(packed,encoding="utf-8",newline="") as text:
            w=csv.writer(text)
            w.writerow(names)
            w.writerows(rows)


def source(temp, out, group, oid, event, tag):
    ingest.immutable(ingest.PLAN)
    raw=ingest.blob(ingest.PLAN,ingest.PREFIX+"/inventory/shards.json")
    need(0<=group<192 and oid in json.loads(raw)[group],"ObsID is outside the original group")
    run=upload_recovery.run_info(ingest.RUN,ingest.ATTEMPT,ingest.CODE)
    item=upload_recovery.get_item(f"chandra-evidence-{group}",run)
    need(item and not item["expired"],"Original producer artifact unavailable")
    file=temp/"producer.zip"
    ingest.api(f"artifacts/{item['id']}/zip",file)
    need("sha256:"+sha(file)==item["digest"] and file.stat().st_size==item["size_in_bytes"],"Changed producer archive")
    with ZipFile(file) as z:
        names=z.namelist()
        need(len(names)==len(set(names)),"Duplicated producer archive member")
        name=f"observations/{oid}.json.gz"
        ledger=z.read(name)
        data=json.loads(gzip.decompress(ledger))
        need(data["obsid"]==oid,"Foreign original observation")
        [p]=[q for q in data["products"] if q["input"]["sha256"]==event]
        [r]=[v for v in data["acquired"] if v["sha256"]==event]
        need(p["input"]==r and p.get("complete") and r["obsid"]==oid,"Incomplete or conflicting product")
        stream,n=tag.split("/")
        need(stream=="real" and str(int(n))==n and int(n)>=0,"Expected a real trigger")
        trig=p["search"]["triggers"][int(n)]
        diag=trig["diagnostics"]
        need(trig["trigger_key"]==tag and trig["obsid"]==oid and diag["input_sha256"]==event,"Foreign trigger")
        base=f"photometry/{oid}/{r['filename']}/{tag}/"
        dbytes=z.read(base+"diagnostics.json")
        need(json.loads(dbytes)==diag,"Producer diagnostic and embedded receipt disagree")
        rows=z.read(base+"event_rows.csv.gz")
        need(digest(rows)==diag["files"]["event_rows.csv.gz"],"Changed archived event-row bytes")
    (out/"source_event_rows.csv.gz").write_bytes(rows)
    save(out/"source.json",{"source_run":ingest.RUN,"source_attempt":ingest.ATTEMPT,"source_code":ingest.CODE,
        "assignment_commit":ingest.PLAN,"assignment_sha256":digest(raw),"group":group,"archive":item,
        "ledger_member":name,"ledger_sha256":digest(ledger),"diagnostics_member":base+"diagnostics.json",
        "diagnostics_sha256":digest(dbytes),"event_rows_member":base+"event_rows.csv.gz","event_rows_sha256":digest(rows),
        "input":r,"original_trigger":trig})
    return r,trig


def verify_rows(evt,trig,expected,temp):
    d,lo,hi=diagnostics._get_data(evt,trig)
    rad=trig["aperture_radius_px"]
    dr=np.hypot(d["x"]-trig["aperture_x"],d["y"]-trig["aperture_y"])
    p0=trig["peak_start_met_s"]
    p1=p0+trig["peak_window_s"]
    rows=temp/"regenerated.csv.gz"
    diagnostics._events_csv(rows,d,dr<=rad,(dr>=2*rad)&(dr<=4*rad),trig["start_met_s"],trig["stop_met_s"],p0,p1)
    need(sha(rows)==sha(expected)==trig["diagnostics"]["files"]["event_rows.csv.gz"],
         "Pinned EVT2 cannot regenerate the deposited row table byte for byte")
    selected={int(x) for x in d["row"][(dr<=rad)&(d["t"]>=p0)&(d["t"]<p1)]}
    need(len(selected)==trig["peak_counts"],"Peak row population differs from original")
    return selected,lo,hi


def status_value(value):
    a=np.asarray(value)
    if a.ndim==0:
        return str(int(a))
    return "".join("1" if b else "0" for b in a.reshape(-1))


def bounds(v):
    return [float(np.min(v)),float(np.max(v))] if len(v) else None


def alignment(cx,cy,x,y):
    def peak(a):
        return max(Counter(int(i) for i in a).values(),default=0)
    return {"chipx":bounds(cx),"chipy":bounds(cy),"sky_x":bounds(x),"sky_y":bounds(y),
            "distinct_chip_pixels":len(set(zip(cx.tolist(),cy.tolist()))),
            "max_on_one_chip_column":peak(cx),"max_on_one_chip_row":peak(cy)}


def analyze(hdu,trig,known,lo,hi,out):
    a=hdu.data
    cols={k.upper():k for k in hdu.columns.names}
    required={"TIME","CCD_ID","CHIPX","CHIPY","X","Y","GRADE","STATUS","ENERGY"}
    need(required<=cols.keys(),"Required level-2 columns missing")
    col=lambda k: np.asarray(a[cols[k]])
    time=col("TIME").astype(float)
    ccd=col("CCD_ID")
    c=int(trig["ccd_id"])
    rows=np.flatnonzero((ccd==c)&(time>=lo)&(time<hi))
    need(len(rows)>0,"Missing source CCD good-time rows")
    t=time[rows]
    x,y,cx,cy=(col(k)[rows].astype(float) for k in ("X","Y","CHIPX","CHIPY"))
    grades=col("GRADE")[rows]
    energy=col("ENERGY")[rows].astype(float)
    statuses=col("STATUS")[rows]
    clean=(statuses==0).reshape(len(rows),-1).all(axis=1)&np.isin(grades,[0,2,3,4,6])&(energy>=500)&(energy<=7000)
    xc,yc=float(trig["aperture_x"]),float(trig["aperture_y"])
    rad=float(trig["aperture_radius_px"])
    need(np.isfinite([rad,xc,yc]).all() and rad>0,"Invalid original aperture")
    dr=np.hypot(x-xc,y-yc)
    inside=dr<=rad
    ann=(dr>=2*rad)&(dr<=4*rad)
    near_radius=max(10*rad,32.)
    near=dr<=near_radius
    p0=trig["peak_start_met_s"]
    pin=(t>=p0)&(t<p0+trig["peak_window_s"])
    peak=clean&inside&pin
    need(set(rows[peak].tolist())==known,"Unfiltered-event map disagrees with all archived selected peak row IDs")
    is_exp="EXPNO" in cols
    group=col("EXPNO")[rows] if is_exp else t
    if is_exp:
        need(np.issubdtype(group.dtype,np.integer) and (group>=0).all(),"Invalid EXPNO; no frame assertion")
    frames,idx=np.unique(group,return_inverse=True)
    nf=len(frames)
    askey=lambda v: int(v) if is_exp else float(v)
    keykind="EXPNO" if is_exp else "identical_TIME_proxy"
    def count(mask):
        return np.bincount(idx[mask],minlength=nf)
    all_counts=count(np.ones(len(rows),bool))
    source=count(inside&clean)
    rawsrc=count(inside)
    peak_counts=count(peak)
    near_counts=count(near)
    tmin=np.full(nf,np.inf)
    tmax=np.full(nf,-np.inf)
    np.minimum.at(tmin,idx,t)
    np.maximum.at(tmax,idx,t)
    ref=np.flatnonzero(near)
    controls=np.full(len(rows),-1,dtype=np.int16)
    control_counts=[]
    centers=[]
    for j in range(8):
        theta=j*np.pi/4
        xj,yj=xc+5*rad*np.cos(theta),yc+5*rad*np.sin(theta)
        sub=ref[np.hypot(x[ref]-xj,y[ref]-yj)<=rad]
        need((controls[sub]<0).all(),"Overlapping control apertures")
        controls[sub]=j
        control_counts.append(np.bincount(idx[sub[clean[sub]]],minlength=nf))
        centers.append([float(xj),float(yj)])
    controls_by_frame=np.stack(control_counts,axis=1)
    peak_frames=set(np.flatnonzero(peak_counts).tolist())
    raw_frame=np.isin(idx,list(peak_frames))
    optional=("NODE_ID","FLTGRADE","PHA","PI","DETX","DETY","TDETX","TDETY")
    opt={k:col(k)[rows] for k in optional if k in cols}
    def event_rows(mask):
        for n in np.flatnonzero(mask):
            yield [int(rows[n]),askey(frames[idx[n]]),float(t[n]),c,float(x[n]),float(y[n]),int(cx[n]),int(cy[n]),
                   float(energy[n]),int(grades[n]),status_value(statuses[n]),int(clean[n]),int(inside[n]),
                   int(ann[n]),int(peak[n]),int(near[n]),int(controls[n]) if controls[n]>=0 else "",
                   *[opt[k][n].item() if k in opt else "" for k in optional]]
    hdr=["evt2_row_zero_based",keykind,"time_met_s","ccd_id","x","y","chipx","chipy","energy_ev","grade",
         "status_as_stored","survey_grade_energy_status_pass","source_aperture","source_annulus","deposited_peak_photon",
         "nearby","control_aperture",*optional]
    csvz(out/"local_level2_rows.csv.gz",hdr,event_rows(near))
    csvz(out/"full_peak_ccd_frames.csv.gz",hdr,event_rows(raw_frame))
    csvz(out/"frame_counts.csv.gz",[keykind,"first_time_met_s","last_time_met_s","whole_ccd_level2","neighborhood_level2",
         "source_level2","source_survey_selected","deposited_peak_count",*[f"control_{j}_survey_selected" for j in range(8)]],
         ([askey(frames[k]),float(tmin[k]),float(tmax[k]),int(all_counts[k]),int(near_counts[k]),int(rawsrc[k]),
           int(source[k]),int(peak_counts[k]),*map(int,controls_by_frame[k])] for k in range(nf)))
    # The S4 serial-readout artifact can occupy one CHIPY row of one node in a single frame.
    # Inspect every row populated by a selected source photon across the whole trial.
    node = (opt["NODE_ID"].astype(int) if "NODE_ID" in opt else (cx.astype(int)-1)//256)
    trial = clean & inside & (t>=trig["start_met_s"]) & (t<trig["stop_met_s"])
    trial_frames=[]
    for k in np.unique(idx[trial]):
        same_frame = idx==k
        src = same_frame & trial
        patterns=[]
        for nd,row_y in sorted(set(zip(node[src].tolist(), cy[src].astype(int).tolist()))):
            row = same_frame & (node==nd) & (cy==row_y)
            inside_trial = row & src
            outside = row & ~inside
            patterns.append({"node":int(nd),"chipy":int(row_y),
                "source_counts":int(inside_trial.sum()),"all_ccd_level2_counts":int(row.sum()),
                "outside_aperture_level2_counts":int(outside.sum()),
                "outside_aperture_survey_selected_counts":int((outside&clean).sum()),
                "source_chipx_range":bounds(cx[inside_trial]),"entire_row_chipx_range":bounds(cx[row]),
                "outside_aperture_chipx_range":bounds(cx[outside])})
        trial_frames.append({"frame_key":askey(frames[k]),"source_trial_counts":int(src.sum()),
                             "source_peak_counts":int(peak_counts[k]),"source_node_rows":patterns})
    frame_detail=[]
    for k in sorted(peak_frames):
        m=peak&(idx==k)
        frame_detail.append({"key":askey(frames[k]),"peak_rows":rows[m].tolist(),"peak_photons":int(m.sum()),
            "whole_ccd_level2_rows":int(all_counts[k]),"local_level2_rows":int(near_counts[k]),
            "time_span_met_s":bounds(t[m]),"peak_geometry":alignment(cx[m],cy[m],x[m],y[m]),
            "grades":grades[m].astype(int).tolist(),"flight_grades":opt["FLTGRADE"][m].astype(int).tolist() if "FLTGRADE" in opt else None,
            "control_counts":controls_by_frame[k].tolist()})
    threshold=int(max(peak_counts,default=0))
    def hist(x):
        return {str(k):v for k,v in sorted(Counter(map(int,x)).items())}
    return {"frame_key":keykind,"observed_ccd_frame_keys":nf,"peak_photons":int(peak.sum()),
            "peak_frame_keys":len(peak_frames),"distinct_peak_time_tags":int(len(np.unique(t[peak]))),
            "max_ccd_intraframe_time_span_s":float(np.max(tmax-tmin)),
            "frame_time_header":{k:hdu.header.get(k) for k in ("TIMEDEL","EXPTIME","TIMEPIXR","TIMESYS","TIMEUNIT")},
            "independent_exposure_numbers_in_peak":len(peak_frames) if is_exp else None,"peak_frames":frame_detail,
            "peak_photon_geometry":alignment(cx[peak],cy[peak],x[peak],y[peak]),"level2_columns":list(cols),
            "absent_optional_columns":sorted(set(optional)-set(opt)),"trigger_gti_met_s":[lo,hi],"instrument_ccd":c,
            "node_identifier":"NODE_ID" if "NODE_ID" in opt else "inferred from 256-column CHIPX intervals",
            "node_row_patterns_across_trial":trial_frames,
            "source_center_sky_px":[xc,yc],"source_radius_px":rad,"local_radius_px":near_radius,
            "control_centers_sky_px":centers,"control_radius_px":rad,"n_control_apertures":8,
            "max_peak_population_per_frame":threshold,"control_frame_aperture_cells":int(controls_by_frame.size),
            "control_cells_at_or_above_peak_frame":int(np.count_nonzero(controls_by_frame>=threshold)) if threshold else None,
            "control_occupancy_histogram":hist(controls_by_frame.ravel()),
            "control_max_occupancy_per_aperture":controls_by_frame.max(axis=0).tolist(),
            "source_occupancy_outside_peak_frames":hist(source[~np.isin(np.arange(nf),list(peak_frames))]),
            "peak_frame_control_histogram":hist(controls_by_frame[sorted(peak_frames)].ravel())}



def brief(m,group,oid,event,tag):
    n=m["peak_photons"]
    k=m["independent_exposure_numbers_in_peak"]
    frame_text=(f"{n} peak photons occupy {k} distinct ACIS exposure number{'s' if k!=1 else ''}"
                if k is not None else f"{n} peak photons occupy {m['distinct_peak_time_tags']} distinct time tags; EXPNO is unavailable")
    geom=m["peak_photon_geometry"]
    lines=[f"# Detector-frame check for Chandra ObsID {oid}", "", f"{frame_text}. "
           f"Their {geom['distinct_chip_pixels']} distinct detector pixels include at most "
           f"{geom['max_on_one_chip_column']} photons on a common CHIPX column and {geom['max_on_one_chip_row']} "
           "on a common CHIPY row.", "",
           f"The brightest source peak frame contains {m['max_peak_population_per_frame']} selected photons. Across "
           f"{m['observed_ccd_frame_keys']} recorded CCD frame keys, "
           f"{m['control_cells_at_or_above_peak_frame']} of {m['control_frame_aperture_cells']} control aperture–frame "
           "pairs reached that count. Eight equal-radius sky apertures lie five aperture radii from the source, "
           "spaced every 45 degrees. Varying exposure, detector edges, diffuse emission, crowding and the choice "
           "to inspect this trigger limit this descriptive comparison. The original survey controls and selection remain the basis for its false-alarm estimate.",
           "", "The event table is pipeline-filtered level 2. Rows in the local and full-frame tables precede this audit’s grade, "
           "energy and status selection. Events rejected by onboard or level-2 processing and unsent pixel-island data are unavailable. "
           "Only CCD frame keys represented by at least one surviving level-2 row enter the denominator. "
           "A shared frame gives no sub-frame arrival ordering. Frame occupancy alone cannot determine particle or celestial origin.",
           "", "[Peak-frame CCD rows](full_peak_ccd_frames.csv.gz) · [Local rows in the full GTI](local_level2_rows.csv.gz) · "
           "[All frame and control counts](frame_counts.csv.gz) · [Measurement and file hashes](result.json) · "
           "[Primary bytes and trigger](source.json) · [Original row table](source_event_rows.csv.gz) · "
           f"[CXC event-column definitions]({DOC})", "",
           f"Primary group {group}, input SHA-256 \x60{event}\x60, trigger \x60{tag}\x60. "
           "The original event-row table was regenerated byte for byte before comparing the frames.", ""]
    if m["instrument_ccd"] == 8:
        lines.extend(["", "## ACIS-S4 serial readout", "",
            "ACIS-S4 can record anomalous events in one frame along one CHIPY row of one readout node. "
            "The row can change between frames. The table includes each row with multiple selected source photons "
            "in one frame anywhere in the trial. The measurement file lists singly occupied rows as well. "
            "Outside counts include every recorded level-2 event in the same node and row beyond the source aperture. "
            f"Node assignments use {m['node_identifier']}. "
            "[CXC S4 streak description](https://cxc.cfa.harvard.edu/ciao/why/destreak.html)", "",
            "| Frame | Node | CHIPY | Source counts | Outside counts | Full CHIPX range |",
            "|---|---:|---:|---:|---:|---|"])
        for frame in m["node_row_patterns_across_trial"]:
            for row in frame["source_node_rows"]:
                if row["source_counts"] >= 2:
                    low,high=row["entire_row_chipx_range"]
                    lines.append(f"| {frame['frame_key']} | {row['node']} | {row['chipy']} | "
                                 f"{row['source_counts']} | {row['outside_aperture_level2_counts']} | "
                                 f"{low:.0f}–{high:.0f} |")
    return "\n".join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--group",type=int,default=GROUP)
    p.add_argument("--obsid",type=int,default=OID)
    p.add_argument("--input-sha256",default=EVENT)
    p.add_argument("--trigger",default=TAG)
    p.add_argument("--retained",type=Path,action="append",default=[])
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO and
         os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH and
         os.getenv("GITHUB_WORKFLOW_REF","").split("@",1)[0]==ingest.REPO+"/"+ingest.FLOW,
         "Use the Chandra research workflow")
    import re
    need(ingest.HEX.fullmatch(a.input_sha256) and re.fullmatch(r"real/(0|[1-9]\d*)",a.trigger) and a.obsid>0,
         "Invalid source or trigger identity")
    output=a.output or campaign.EVIDENCE/"frame_forensics"/str(a.obsid)/a.input_sha256/a.trigger
    # Preserve earlier attempts and outputs from any later changed implementation.
    out=output/f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    out.mkdir(parents=True,exist_ok=True)
    need(not (out/"result.json").exists(),"A result for this attempt already exists")
    attempts=[]
    try:
        with tempfile.TemporaryDirectory(dir=os.getenv("RUNNER_TEMP")) as place:
            temp=Path(place)
            rec,trig=source(temp,out,a.group,a.obsid,a.input_sha256,a.trigger)
            file=followup.event_bytes(rec,a.obsid,temp,a.retained,attempts)
            evt=Path(file["path"])
            known,lo,hi=verify_rows(evt,trig,out/"source_event_rows.csv.gz",temp)
            if (a.group,a.obsid,a.input_sha256,a.trigger)==(GROUP,OID,EVENT,TAG):
                need(len(known)==6,"The original target does not contain its six pinned peak rows")
            with fits.open(evt,memmap=False) as hdus:
                m=analyze(hdus["EVENTS"],trig,known,lo,hi,out)
        doc=brief(m,a.group,a.obsid,a.input_sha256,a.trigger)
        (out/"README.md").write_text(doc)
        save(out/"acquisition.json",attempts)
        files=("source.json","source_event_rows.csv.gz","frame_counts.csv.gz","local_level2_rows.csv.gz",
               "full_peak_ccd_frames.csv.gz","acquisition.json","README.md")
        save(out/"result.json",{"status":"complete","target":{"group":a.group,"obsid":a.obsid,"input_sha256":a.input_sha256,
             "trigger":a.trigger},"source_run":ingest.RUN,"execution_run":os.environ["GITHUB_RUN_ID"],
             "execution_attempt":os.environ["GITHUB_RUN_ATTEMPT"],"execution_commit":os.environ["GITHUB_SHA"],
             "implementation_sha256":sha(Path(__file__)),"dependencies_sha256":{k:sha(campaign.HERE/k) for k in ("diagnostics.py","archive.py","followup.py","ingest_primary.py","requirements.txt")},"regenerated_rows_sha256":sha(out/"source_event_rows.csv.gz"),
             "files":{k:sha(out/k) for k in files},"measurement":m})
        print(doc.splitlines()[2],flush=True)
    except Exception as exc:
        save(out/"acquisition.json",attempts)
        save(out/"failure.json",{"status":"incomplete","type":type(exc).__name__,"message":str(exc),
             "identity":[a.group,a.obsid,a.input_sha256,a.trigger],"implementation_sha256":sha(Path(__file__))})
        raise


if __name__=="__main__":
    main()
