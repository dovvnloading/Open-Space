#!/usr/bin/env python3
"""A single bounded checkpoint of full HATS publisher source rows.
Original strict-selected JSONLs CANNOT recover rows excluded by old cuts.
Reuse original 768 S3 key/ETag receipts, stream projected Parquet row groups,
store exact global eligibility counts and rank-preserving top 200/partition.
"""
import gzip,hashlib,json,math,sys,os
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import pyarrow.parquet as pq
import pyarrow.fs as pafs
import selection

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
OLD=ROOT.parent/"euclid-dr1-wide-highz-vetting"
sys.path.insert(0,str(OLD))
import scan
BATCH=24

def sha(b):return hashlib.sha256(b).hexdigest()
def js(path,x):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+"\n")
def finite(x):
    try:
        f=float(x)
        return f if math.isfinite(f) else None
    except (TypeError,ValueError):return None

def score_key(r):return (-r["rank_score"],r["id"])
def compact_row(tab,i,arrays,part):
    # Store enough for future optical/flag/motion assessment; no hidden cut.
    v=tab.slice(int(i),1).to_pylist()[0]
    def scalar(k):return scan.scalar(v.get(k))
    def number(k):return finite(arrays[k][i])
    return {
        "id":str(scalar("object_id")),"ra":scalar("ra"),"dec":scalar("dec"),
        "tileid":scalar("tileid"),"Norder":part["Norder"],"Npix":part["Npix"],
        "rank_score":round(number("score"),6),
        "Y_minus_J_AB":number("YJ"),"Y_minus_J_2sig_limit":number("YJ_LIMIT"),
        "Y_minus_J_used":number("YJ_USED"),"J_minus_H_AB":number("JH"),
        "VIS_SNR":round(number("VSNR"),4),
        "J_SNR":round(number("JSNR"),4),
        "H_SNR":round(number("HSNR"),4),
        "optical_max_SNR":number("optical_max"),
        "ground_optical_gt3":bool(arrays["optical_gt3"][i]),
        "nonzero_J_H_flags":bool(arrays["bad_flags"][i]),
        "vis_det":scalar("vis_det"),"flag_vis":scalar("flag_vis"),
        "flag_j":scalar("flag_j"),"flag_h":scalar("flag_h"),
        "det_quality_flag":scalar("det_quality_flag"),
        "mumax_minus_mag":scalar("mumax_minus_mag"),
        "ellipticity":scalar("ellipticity"),
        "spurious_prob":scalar("spurious_prob"),
        "source_s3_key":part["key"],"source_etag":part["etag"],
        "image_artifacts_cleared":False,"spectroscopically_new_Q":False
    }

def scan_one(part,fs):
    # Immutable source receipt from fully replayed PR208 output, not re-listing
    # against a possibly changing HATS prefix.
    key=part["key"]
    head=scan.SESSION.head(scan.BASE+key,timeout=60)
    head.raise_for_status()
    if head.headers.get("ETag")!=part["etag"]:
        raise RuntimeError("Source ETag changed: "+key)
    with fs.open_input_file(scan.BUCKET+"/"+key) as stream:
        parquet=pq.ParquetFile(stream)
        if parquet.metadata.num_rows!=part["rows"]:
            raise RuntimeError("Parquet rows changed: "+key)
        cols=[x for x in scan.COLUMNS if x in parquet.schema_arrow.names]
        missing=scan.REQUIRED-set(cols)
        if missing or "mumax_minus_mag" not in cols:
            raise RuntimeError("Missing required data fields "+str(missing))
        selected=0; clean_opt=0; flagged=0; total=0; top=[]
        for rg in range(parquet.num_row_groups):
            tab=parquet.read_row_group(rg,columns=cols,use_threads=False)
            a=selection.arrays(tab)
            mask=a["eligible"]
            ix=np.flatnonzero(mask)
            total+=len(mask)
            selected+=len(ix)
            clean_opt+=int(np.count_nonzero(mask & ~a["optical_gt3"]))
            flagged+=int(np.count_nonzero(mask & a["bad_flags"]))
            if not len(ix):continue
            # Each row group retains enough for global top-200 correctness.
            if len(ix)>selection.TOP_PER_PARTITION:
                q=selection.TOP_PER_PARTITION
                ix=ix[np.argpartition(a["score"][ix],-q)[-q:]]
            top.extend(compact_row(tab,int(i),a,part) for i in ix)
            top.sort(key=score_key)
            top=top[:selection.TOP_PER_PARTITION]
        if total!=part["rows"]:
            raise RuntimeError("Projected row group count mismatch "+key)
        return {"Norder":part["Norder"],"Npix":part["Npix"],
                "key":key,"etag":part["etag"],"read_rows":total,
                "eligible_rows":selected,"eligible_optical_not_gt3":clean_opt,
                "eligible_flagged_JH":flagged,
                "retained_top_rows":len(top),
                "rank_truncation_threshold":top[-1]["rank_score"] if len(top)==selection.TOP_PER_PARTITION else None,
                "top":top}

def prior_partitions():
    seq=[]
    for batch in range(32):
        p=OLD/f"output/batches/batch_{batch:03d}_manifest.json"
        m=json.loads(p.read_text())
        if m["index_sha256"]!=scan.FROZEN_INDEX_SHA256:
            raise RuntimeError("Original source index SHA mismatch")
        seq.extend(m["partitions"])
    if len(seq)!=768 or sum(r["rows"] for r in seq)!=scan.PUBLISHED_ROWS:
        raise AssertionError("Original exact HATS source coverage changed")
    return seq

def main():
    if not (OUT/"calibration.json").exists():raise RuntimeError("Control replay first")
    cal=json.loads((OUT/"calibration.json").read_text())
    if cal["matched_quasar_recovery"]["numerator"]!=21 or cal["matched_dwarf_retention"]["numerator"]!=4:
        raise RuntimeError("Unfrozen selector")
    (OUT/"batches").mkdir(parents=True,exist_ok=True)
    parts=prior_partitions()
    completed=[]
    for b in range(32):
        manifest=OUT/f"batches/batch_{b:03d}_manifest.json"
        candidates=OUT/f"batches/batch_{b:03d}_top.jsonl.gz"
        if manifest.exists():
            m=json.loads(manifest.read_text())
            if not candidates.exists() or sha(candidates.read_bytes())!=m["top_sha256"]:
                raise RuntimeError("Checkpoint digest mismatch batch "+str(b))
            completed.append(b)
    if completed != list(range(len(completed))):
        raise RuntimeError("Missing batch before later completed batch")
    if len(completed)==32:
        print("RESELECTION_COMPLETE_768_PARTITIONS",flush=True)
        return
    b=len(completed)
    sub=parts[b*BATCH:(b+1)*BATCH]
    fs=pafs.S3FileSystem(anonymous=True,region="us-east-1",request_timeout=120,connect_timeout=25)
    found={}
    with ThreadPoolExecutor(max_workers=int(os.getenv("EUCLID_PARTITION_WORKERS","4"))) as pool:
        futures={pool.submit(scan_one,p,fs):i for i,p in enumerate(sub)}
        for f in as_completed(futures):
            i=futures[f]
            x=f.result()
            found[i]=x
            print("PARTITION",b,i+1,"rows",x["read_rows"],"new_eligible",x["eligible_rows"],flush=True)
    if len(found)!=len(sub):raise RuntimeError("Incomplete partition batch")
    ordered=[found[i] for i in range(len(sub))]
    rows=[row for part in ordered for row in part.pop("top")]
    # Stable serialized gzip; hash and register this checkpoint.
    cp=OUT/f"batches/batch_{b:03d}_top.jsonl.gz"
    with gzip.open(cp,"wt",encoding="utf8",compresslevel=6) as f:
        for r in rows:f.write(json.dumps(r,sort_keys=True,allow_nan=False)+"\n")
    m={"status":"COMPLETE","batch":b,
       "selection_id":json.loads((ROOT/"freeze.json").read_text())["selection_id"],
       "index_sha256":scan.FROZEN_INDEX_SHA256,
       "partitions":ordered,
       "read_rows":sum(p["read_rows"] for p in ordered),
       "eligible_rows":sum(p["eligible_rows"] for p in ordered),
       "retained_top_rows":len(rows),
       "top_sha256":sha(cp.read_bytes())}
    js(OUT/f"batches/batch_{b:03d}_manifest.json",m)
    js(OUT/"status.json",{
       "phase":"RESELECTING_PUBLISHER_WIDE_HATS",
       "completed_batches":b+1,"total_batches":32,
       "scanned_partitions":(b+1)*BATCH,
       "selected_count_in_completed_batches":sum(
           json.loads((OUT/f"batches/batch_{j:03d}_manifest.json").read_text())["eligible_rows"]
           for j in range(b+1)),
       "calibrated_full_footprint_completeness":False})
    print("BATCH_FINISHED",b+1,"/32 total",m["eligible_rows"],flush=True)
if __name__=="__main__":main()
