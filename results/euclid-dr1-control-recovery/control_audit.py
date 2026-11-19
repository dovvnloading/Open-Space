#!/usr/bin/env python3
"""Exact per-object failure audit for catalogue-matched Euclid z>6.5 controls.
Only networked GitHub Actions runs this module. Frozen source ETags are checked.
"""
import hashlib
import json
import math
import sys
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pyarrow.parquet as pq
import pyarrow.fs as pafs

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "euclid-dr1-wide-highz-vetting"
OUT = HERE / "output"
sys.path.insert(0, str(OLD))
import scan

def save(name, x):
    OUT.mkdir(exist_ok=True)
    (OUT / name).write_text(json.dumps(x, indent=2, sort_keys=True, allow_nan=False) + "\n")

def val(x):
    try:
        x = float(x)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None

def ratio(a, b):
    a, b = val(a), val(b)
    return a / b if a is not None and b is not None and b > 0 else None

def reasons(r):
    """Independently replay the old selection's conjuncts, not broad alone."""
    get = lambda k: val(r.get(k))
    v, ve = get("flux_vis_3fwhm_aper"), get("fluxerr_vis_3fwhm_aper")
    y, ye = get("flux_y_3fwhm_aper"), get("fluxerr_y_3fwhm_aper")
    j, je = get("flux_j_3fwhm_aper"), get("fluxerr_j_3fwhm_aper")
    h, he = get("flux_h_3fwhm_aper"), get("fluxerr_h_3fwhm_aper")
    ell = get("ellipticity")
    spur = get("spurious_prob")
    vd = get("vis_det")
    fj, fh = get("flag_j"), get("flag_h")
    snv, snj, snh = ratio(v,ve), ratio(j,je), ratio(h,he)
    yj = (-2.5 * math.log10(y/j) if y is not None and j is not None and y>0 and j>0 else None)
    jh = (-2.5 * math.log10(j/h) if j is not None and h is not None and j>0 and h>0 else None)
    ylower = (-2.5 * math.log10(max(2*ye,y if y is not None else -math.inf)/j)
              if ye is not None and j is not None and j>0 and ye>=0 and max(2*ye,y if y is not None else -math.inf)>0 else None)
    tests = {
        "VIS_DET_equals_zero": vd == 0,
        "J_flux_gt_0p65": j is not None and j > .65,
        "H_flux_positive": h is not None and h > 0,
        "Y_less_than_0p85J": y is not None and j is not None and y < .85*j,
        "J_greater_than_0p65H": j is not None and h is not None and j > .65*h,
        "J_error_positive": je is not None and je > 0,
        "H_error_positive": he is not None and he > 0,
        "spurious_lt_0p7_or_missing": spur is None or spur < .7,
        "J_SNR_ge_6": snj is not None and snj >= 6,
        "H_SNR_ge_4": snh is not None and snh >= 4,
        "VIS_error_positive": ve is not None and ve > 0,
        "VIS_abs_SNR_lt_2": snv is not None and abs(snv) < 2,
        "YJ_ge_0p4_or_limit": (yj is not None and yj >= .4) if y is not None and y > 0 else
                (ylower is not None and ylower >= .4),
        "JH_le_0p45": jh is not None and jh <= .45,
        "ellipticity_lt_0p30": ell is not None and ell < .30,
        "J_flag_zero": fj == 0,
        "H_flag_zero": fh == 0,
    }
    derived = {"VIS_SNR_3fwhm":snv, "J_SNR_3fwhm":snj, "H_SNR_3fwhm":snh,
               "Y_minus_J_AB":yj, "J_minus_H_AB":jh, "Y_minus_J_2sig_limit":ylower}
    return tests, derived

def get_tasks():
    controls = json.loads((OLD / "output/controls.json").read_text())
    selected = {c["name"]: c for c in controls if c["nearest"]}
    locations = {}
    for batch in range(32):
        manifest = json.loads((OLD / f"output/batches/batch_{batch:03d}_manifest.json").read_text())
        for part in manifest["partitions"]:
            for hit in part["controls"]:
                name = hit["name"]
                if name not in selected: continue
                if str(hit["object_id"]) != str(selected[name]["nearest"]["object_id"]): continue
                identity = (part["Norder"], part["Npix"], part["key"], part["etag"])
                if name in locations and locations[name][0] != identity:
                    raise AssertionError("Same nearest object in different source partition: "+name)
                locations[name] = (identity, hit)
    missing = sorted(set(selected) - set(locations))
    if missing: raise AssertionError("Frozen nearest controls lack manifest source: " + repr(missing))
    groups = defaultdict(list)
    for name,(identity,hit) in locations.items():
        groups[identity].append((name, selected[name], hit))
    return controls, groups

def audit_partition(identity, items):
    order,pix,key,expected_etag=identity
    head = scan.SESSION.head(scan.BASE + key, timeout=60)
    head.raise_for_status()
    actual = head.headers.get("ETag")
    if actual != expected_etag:
        raise RuntimeError("S3 product ETag changed: "+key)
    fs = pafs.S3FileSystem(anonymous=True, region="us-east-1",
                          request_timeout=100, connect_timeout=20)
    targets=defaultdict(list)
    for name,c,hit in items:
        targets[int(hit["object_id"])].append((name,c,hit))
    found={}
    with fs.open_input_file(scan.BUCKET+"/"+key) as fh:
        parquet=pq.ParquetFile(fh)
        cols=[k for k in scan.COLUMNS if k in parquet.schema_arrow.names]
        for rg in range(parquet.num_row_groups):
            ids=parquet.read_row_group(rg,columns=["object_id"]).column(0).to_numpy(zero_copy_only=False)
            indices=np.flatnonzero(np.isin(ids, list(targets)))
            if not len(indices): continue
            table=parquet.read_row_group(rg,columns=cols)
            for idx in indices:
                source=table.slice(int(idx),1).to_pylist()[0]
                source={k:scan.scalar(v) for k,v in source.items()}
                for name,c,hit in targets[int(ids[idx])]:
                    checks,derived=reasons(source)
                    if name in found: raise AssertionError("Control object occurs twice: "+name)
                    found[name]={
                        "name":name,"class":c["class"],"z":c["z"],"source_catalog":c["source"],
                        "object_id":str(hit["object_id"]),"sep_arcsec":hit["sep_arcsec"],
                        "original_strict_pass":hit["strict_pass"],"original_broad_pass":hit["broad_pass"],
                        "source_s3_key":key,"source_etag":actual,
                        "Norder":order,"Npix":pix,"row_group":rg,
                        "source":source,"derived":derived,"original_cut_checks":checks,
                        "failed_original_cuts":[k for k,v in checks.items() if not v],
                        "independent_strict_pass":all(checks.values())
                    }
    wanted={n for n,_,_ in items}
    if set(found)!=wanted: raise AssertionError("Missing object IDs in "+key+": "+repr(wanted-set(found)))
    for name,c,hit in items:
        if found[name]["independent_strict_pass"]!=bool(hit["strict_pass"]):
            raise AssertionError("Independent replay mismatch "+name)
    return found

def main():
    controls,groups=get_tasks()
    print("SOURCE_PARTITIONS",len(groups),"MATCHED_CONTROLS",sum(len(t) for t in groups.values()),flush=True)
    output={}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(audit_partition,ident,items):ident for ident,items in groups.items()}
        for future in as_completed(futures):
            output.update(future.result())
            print("AUDITED",len(output),"/",sum(map(len,groups.values())),flush=True)
    rows=[]
    for c in controls:
        v=output.get(c["name"])
        rows.append(v if v is not None else
                   {"name":c["name"],"class":c["class"],"z":c["z"],
                    "status":"NO_3ARCSEC_CATALOG_MATCH_COVERAGE_UNKNOWN",
                    "source_catalog":c["source"]})
    matched_q=[r for r in rows if r["class"]=="quasar" and "source" in r]
    recovered=sum(r["independent_strict_pass"] for r in matched_q)
    dwarfs=[r for r in rows if r["class"]=="dwarf" and "source" in r]
    counts=Counter(key for r in matched_q if not r["independent_strict_pass"]
                   for key in r["failed_original_cuts"])
    summary={"status":"CONTROL_RAW_REPLAY_COMPLETE",
             "matched_Q":len(matched_q),"old_strict_recovered_Q":recovered,
             "matched_dwarfs":len(dwarfs),"matched_control_sources":len(output),
             "source_partitions":len(groups),
             "failed_cut_frequencies_nonexclusive":dict(counts),
             "unmatched_Q_denominator_excluded":sum(r["class"]=="quasar" and "source" not in r for r in rows),
             "footprint_completeness":None, "discovery_level":False}
    save("control_audit.json",rows)
    save("audit_summary.json",summary)
    print("AUDIT_SUMMARY",json.dumps(summary),flush=True)

if __name__=="__main__":
    main()
