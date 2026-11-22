"""Rebuild the imaging-scout wavelength shuffles on their original assigned positions."""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import gzip
import json

import numpy as np

from acquisition import Spectra, prepare
from desi_lines import external_lines
from evidence_summary import ledger
from ew_imaging_followup import decision, digest, support, trim, write
from legacy_bricks import legacy_brick_cones
import spectral as spec
from vetting import Ancillary, Cache, atom, fixed_external_pairs, PAIR_EXTERNAL, PAIR_REPEAT, validate_local


def reconstruct(cat, register, counts):
    """The frozen register records the sorted target array used for each permutation."""
    indices=np.array([row["row_index"] for _,row in register],dtype=int)
    if np.any(np.diff(indices)<=0):
        raise ValueError("Control register is not in original strictly increasing target order")
    for first in range(0,len(indices),64):
        block=register[first:first+64]
        local=[j for j,(_,r) in enumerate(block) if r["status"]=="science_peak_for_imaging"]
        if not local:
            continue
        ix=indices[first:first+64]
        f=np.array(cat.flux[ix],dtype=np.float32,order="C")
        e=np.array(cat.error[ix],dtype=np.float32,order="C")
        e[(f==-999)|~np.isfinite(f)|~np.isfinite(e)|(e<=0)]=np.nan
        seed=572025+first
        y,err=spec.shuffle_fibres(f,e,seed=seed)
        recovered=defaultdict(list)
        for j,p in spec.search_batch(cat.wave,y[local],err[local],density_unit=1e-17):
            recovered[local[j]].append(p)
        for j in local:
            line,row=block[j]
            tests=[p for p in recovered[j] if abs(p["wave"]-row["wave"])<1e-5]
            keys=("wave","flux","flux_err","snr","sigma_A","noise_scale")
            if len(tests)!=1 or any(not np.isclose(tests[0][k],row[k],rtol=1e-6,atol=0) for k in keys):
                raise ValueError(f"Checksum-verified source files could not reproduce scout control at ledger line {line}")
            counts["verified_frozen_primary_fits"]+=1
            m=cat.metadata(int(ix[j]))
            if any(row[k]!=m[k] for k in ("shotid","detectid","source_id")):
                raise ValueError("Catalog row identity differs from frozen control")
            ra=m.get("ra_det") if m.get("ra_det") is not None else row["ra"]
            dec=m.get("dec_det") if m.get("dec_det") is not None else row["dec"]
            p=dict(row,ra=ra,dec=dec,channel="ew",stream="imaging_rescue",release=cat.release,
                   shuffled_detail_complete=True,row_is_synthetic=True)
            p.update(spec.detail(cat.wave,y[j],err[j],p,density_unit=1e-17))
            design=dict(seed=seed,block_offset=first,recipient_offset=j,donor_rows=ix.tolist())
            yield line,{**m,**p,"catalog":m},design


def archive_at_positions(out,positions):
    cones={}
    bands=defaultdict(list)
    if not positions:
        return cones,bands
    for n,record in ledger(out/"native.jsonl.gz"):
        cone=record["cone"]
        key=(cone.get("ra"),cone.get("dec"))
        if key not in positions or cone["status"] not in ("queried","outside_footprint"):
            continue
        if key not in cones or cones[key][1]["status"]!="queried":
            cones[key]=(n,cone)
        if cone["status"]=="queried":
            bands[(key,cone.get("selected_release"),cone.get("selected_cap"))].append((n,record["bands"]))
    return cones,bands


def run_control(indir,out,release,cache_root):
    dest=out/"scout-control"
    dest.mkdir(parents=True,exist_ok=True)
    source=indir/"ew-rescue-control-spectra.jsonl.gz"
    manifest=indir/"ew-rescue-control-summary.json"
    registers=list(ledger(source))
    old=json.loads(manifest.read_text())["counts"]
    cases=Counter(row["status"] for _,row in registers)
    counts=Counter(eligible_targets=old["eligible_source_observations"],shuffled_spectra=len(registers),
                   singleton_omitted=old.get("control_singleton_not_shuffled",0),
                   preliminary_science_peaks=cases["science_peak_for_imaging"])
    if (counts["eligible_targets"]!=len(registers)+counts["singleton_omitted"] or
        len(registers)!=old["spectra_measured"] or
        any(n!=old.get(k,0) for k,n in cases.items())):
        raise ValueError("Control register and frozen scout summary disagree")
    previous={}
    for n,entry in ledger(indir/"ew-rescue-control-imaging.jsonl.gz"):
        p=entry["source"]
        previous[p["row_index"],p["wave"],p["ra"],p["dec"]]=(n,entry)
    path=prepare(release,cache_root,dest)
    cache=Cache(cache_root)
    aux=Ancillary(cache)
    try:
        with Spectra(path,release) as cat, gzip.open(dest/"assessments.jsonl.gz","wt",compresslevel=3) as audit, \
             gzip.open(dest/"native.jsonl.gz","wt",compresslevel=3) as native_log:
            rows=[]
            for n,rec,design in reconstruct(cat,registers,counts):
                v=validate_local(rec,cat,aux)
                if v["disposition"]=="unresolved":
                    if any(r.get("snr",0)>=4 for r in v.get("secondary_repeats",[])):
                        v["passed"].append(PAIR_REPEAT)
                    if fixed_external_pairs(rec,v):
                        v["passed"].append(PAIR_EXTERNAL)
                rows.append((n,rec,design,v))
            def extra(item):
                _,rec,_,v=item
                if v["disposition"]=="unresolved" and not support(rec,v):
                    v["SPARCL"]=external_lines(rec["ra"],rec["dec"],rec["wave"],"ew",cache.root)
                    if fixed_external_pairs(rec,v) and PAIR_EXTERNAL not in v["passed"]:
                        v["passed"].append(PAIR_EXTERNAL)
            with ThreadPoolExecutor(max_workers=3) as pool:
                list(pool.map(extra,rows))
            positions={(r["ra"],r["dec"]) for _,r,_,v in rows if v["disposition"]=="unresolved" and support(r,v)}
            cones,bands=archive_at_positions(out,positions)
            native_line=0
            rows.sort(key=lambda row:(round(row[1]["ra"],1),round(row[1]["dec"],1),row[1]["ra"]))
            for n,r,design,v in rows:
                key=(r["ra"],r["dec"])
                assessment=dict(source_ledger=source.name,source_line=n,design=design,
                                peak={k:r.get(k) for k in ("source_id","shotid","detectid","row_index","release","ra","dec",
                                      "wave","wave_vac","z_lya","flux","flux_err","snr","sigma_A","ew0","ew0_lower",
                                      "noise_scale","rescue_population","row_is_synthetic","shuffled_detail_complete",
                                      "low_z","lya_family")},local_disposition=v["disposition"])
                if v.get("SPARCL"):
                    counts["additional_external_cone_"+v["SPARCL"].get("status","unknown")]+=1
                if v["disposition"]=="unresolved" and support(r,v):
                    counts["independently_supported"]+=1
                    prior=previous.get((r["row_index"],r["wave"],*key))
                    if prior:
                        ln,entry=prior
                        v.update({k:d for k,d in entry["evidence"].items() if k in ("HSC","Legacy_forced","cube") and d})
                        v["imaging_EW"]=entry["imaging"]
                        assessment["previous_mock_imaging_line"]=ln
                    ref,cone=cones.get(key,(None,None))
                    if cone is None:
                        cone=legacy_brick_cones([key],cache)[0]
                        cones[key]=(None,cone)
                    elif ref is not None:
                        counts["archived_cones_reused"]+=1
                    available=bands[key,cone.get("selected_release"),cone.get("selected_cap")]
                    v,measurements,img=decision(r,v,cone,cache,cat,aux,counts,available)
                    native_line+=1
                    write(native_log,dict(row_index=r["row_index"],wave=r["wave"],source_line=n,
                                          prior_native_cone_line=ref,cone=cone,bands=measurements,imaging=img))
                    assessment["native_line"]=native_line
                    assessment["independently_supported"]=True
                    trim(cache)
                else:
                    assessment["independently_supported"]=False
                    if v["disposition"]=="unresolved":
                        v["reasons"].append("independent redshift or line corroboration absent")
                counts["disposition_"+v["disposition"]]+=1
                write(audit,dict(assessment,vet=v))
    finally:
        aux.close()
    if counts["preliminary_science_peaks"]!=counts["verified_frozen_primary_fits"]:
        raise ValueError("Incomplete frozen scout replay")
    summary=dict(release=release,counts=counts,first_pass_status_counts=cases,catalog_errors=aux.errors,
                 ancillary_and_cube_receipts=cache.inputs,
                 sha256={str(p.relative_to(indir)):digest(p) for p in
                         (source,manifest,indir/"ew-rescue-control-imaging.jsonl.gz",out/"native.jsonl.gz")},
                 method="Position-assigned mocks use the original paired flux/error wavelength permutation across sorted "
                        "64-target groups, seed 572025 plus group offset. Every retained primary flux and error is verified "
                        "against the frozen blind fit. Assigned sky positions keep the real DESI or separate-shot confirmation, "
                        "cube and imaging tests. Earlier photometry reuses the same control or a native continuum measurement "
                        "at the identical position and cube PSF. No observed source-line fit substitutes for a mock primary. "
                        "An additional SPARCL cone is attempted for each locally unresolved control lacking corroboration.")
    (dest/"summary.json").write_text(json.dumps(atom(summary),indent=2,allow_nan=False)+"\n")
    print(release,"dedicated imaging-scout null:",dict(counts),flush=True)
