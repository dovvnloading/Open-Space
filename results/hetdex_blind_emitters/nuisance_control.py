"""Same-shot controls preserving per-source continua, empirical noise and short covariance."""
from collections import Counter
from itertools import groupby
import gzip
import json
from pathlib import Path
import random

import numpy as np

import spectral as spec
from repeats import spectra_for_rows

SEED = 3939941
BLOCK = 4  # native 2-A channels per contiguous donor block
DONORS = 8
STRATUM = 64


def _model(y,e):
    """Piecewise-interpolated robust continuum and whitened residuals."""
    ok = np.isfinite(y) & np.isfinite(e) & (e>0)
    x = np.arange(len(y))
    cuts = np.arange(0,len(y),32)
    center = []
    value = []
    for first in cuts:
        use = ok[first:first+32]
        if use.sum() >= 12:
            center.append(np.mean(x[first:first+32][use]))
            value.append(np.median(y[first:first+32][use]))
    if len(value)<4 or ok.sum()<192:
        return None
    cont = np.interp(x,center,value)
    std = np.divide(y-cont,e,out=np.full(len(y),np.nan),where=ok)
    center_res = np.median(std[ok])
    gain = 1.4826*np.median(np.abs(std[ok]-center_res))
    if not np.isfinite(gain) or gain <= .05:
        return None
    scale = np.median(e[ok])
    colors=[]
    for window in np.array_split(x,3):
        good=window[ok[window]]
        if len(good)<24:
            return None
        colors.append(float(np.arcsinh(np.median(cont[good])/scale)))
    descriptor=np.array([*colors,np.log(scale),np.log(gain),np.mean(ok)])
    return cont,(std-center_res)/gain,gain,descriptor,ok


def _draw(target, donors, err, rngseed):
    """Repeated calls with these inputs reconstruct the identical null."""
    base,_,gain,_,mask = target
    z = np.stack([item[1] for item in donors])
    valid = np.isfinite(z)
    n=len(mask)
    rng=np.random.default_rng(rngseed)
    phase=int(rng.integers(BLOCK))
    bounds=np.r_[0,np.arange(BLOCK-phase,n,BLOCK),n]
    # Random order independently per block. Among donors with the greatest
    # usable overlap choose uniformly, keeping one identity over the block.
    block_index=np.searchsorted(bounds[1:],np.arange(n),side="right")
    cover=np.add.reduceat((valid & mask[None,:]).astype(np.int16),bounds[:-1],axis=1)
    chosen=np.argmax(cover+rng.random(cover.shape),axis=0)
    indices=chosen[block_index]
    sample=z[indices,np.arange(n)]
    good=mask & np.isfinite(sample)
    f=np.full(n,np.nan)
    f[good]=base[good]+gain*err[good]*sample[good]
    return f,np.where(good,err,np.nan)


def _signature(row, pool, descriptors, xyz, ids):
    distance=np.linalg.norm(xyz[pool]-xyz[row],axis=1)
    possible=np.flatnonzero((ids[pool]!=ids[row]) & (distance>np.deg2rad(10/3600)))
    target=descriptors[row]
    scored=[]
    for q in possible:
        donor=int(pool[q])
        desc=descriptors.get(donor)
        if desc is None:
            continue
        d=np.abs(desc-target)
        # Reject extreme mismatches even if the nearest donor in a sparse
        # shot is distant: factor 3 in error/gain and comparable broad color.
        if np.max(d[:3])>1.5 or max(d[3:5])>1.1 or d[5]>.25:
            continue
        score=float(np.sqrt(np.mean(d[:3]**2)+d[3]**2+.5*d[4]**2+d[5]**2))
        scored.append((score,donor))
    scored.sort()
    return scored[:DONORS]


def _json(obj):
    from scan import dump
    return dump(obj)


def _shots(cat):
    order=np.argsort(cat.shot_ids,kind="stable")
    ids=cat.shot_ids[order]
    cuts=np.r_[0,np.flatnonzero(np.diff(ids))+1,len(order)]
    for p,q in zip(cuts[:-1],cuts[1:]):
        if ids[p] > 0:
            yield int(ids[p]),order[p:q]


def _quality_summary(real,fake):
    def measures(obj):
        v=obj[np.isfinite(obj)]
        return float(np.median(v)) if len(v) else None
    return {k:{"real_median":measures(np.asarray(real[k],float)),
               "mock_median":measures(np.asarray(fake[k],float))} for k in real}


def _quality(y,e,model):
    if model is None:
        return (np.nan,np.nan)
    z=model[1]
    ok=np.isfinite(z[1:]) & np.isfinite(z[:-1])
    if ok.sum()<16:
        return (model[2],np.nan)
    # Robust winsorization keeps one emission line from dominating covariance.
    a=np.clip(z[1:][ok],-4,4)
    b=np.clip(z[:-1][ok],-4,4)
    return (model[2],float(np.corrcoef(a,b)[0,1]))


def run_control(cat,out):
    from scan import screen_chunk
    from repeats import vectors
    out=Path(out)
    xyz=vectors(np.asarray(cat.info["ra"],float),np.asarray(cat.info["dec"],float))
    totals=Counter()
    observed=Counter()
    mocked=Counter()
    real_quality={"gain":[],"lag1_correlation":[]}
    mock_quality={"gain":[],"lag1_correlation":[]}
    with gzip.open(out/"nuisance-design.jsonl.gz","wt",compresslevel=3) as design, \
         gzip.open(out/"nuisance-leads.jsonl.gz","wt",compresslevel=3) as leads, \
         gzip.open(out/"nuisance-real-leads.jsonl.gz","wt",compresslevel=3) as real_leads:
        for shot,members in _shots(cat):
            seed=SEED+(1 if cat.release=="sc2" else 0)
            rng=np.random.default_rng([seed,shot])
            nr, remainder=divmod(len(members),STRATUM)
            nr+=int(rng.random()<remainder/STRATUM)
            totals["shot_strata"]+=1
            if not nr:
                continue
            targets=rng.choice(members,size=nr,replace=False)
            pool=rng.choice(members,size=min(len(members),max(32,min(128,5*nr))),replace=False)
            all_rows=np.unique(np.r_[pool,targets])
            f,e=spectra_for_rows(cat,all_rows)
            indexed={int(row):j for j,row in enumerate(all_rows)}
            models={int(row):_model(f[j],e[j]) for j,row in enumerate(all_rows)}
            descriptors={k:m[3] for k,m in models.items() if m is not None}
            plans=[]
            fy,fe,oy,oe=[],[],[],[]
            for target0 in targets:
                target=int(target0)
                totals["drawn"]+=1
                mm=models[target]
                if mm is None:
                    totals["recipient_unusable"]+=1
                    continue
                scores=_signature(target,pool,descriptors,xyz,cat.src_ids)
                if len(scores)<3:
                    totals["fewer_than_three_comparable_donors"]+=1
                    continue
                donorrows=[j for _,j in scores]
                rngseed=seed+1009*target
                t=indexed[target]
                fake,err=_draw(mm,[models[j] for j in donorrows],e[t],rngseed)
                if np.count_nonzero(np.isfinite(err))<192:
                    totals["insufficient_common_bins"]+=1
                    continue
                plan={"release":cat.release,"shotid":shot,"row_index":target,"donor_rows":donorrows,
                      "draw_seed":rngseed,"worst_donor_distance":scores[-1][0],
                      "recipient_gain":float(mm[2]),"usable_recipient_bins":int(mm[4].sum()),
                      "usable_mock_bins":int(np.count_nonzero(np.isfinite(err)))}
                design.write(_json(plan)+"\n")
                fy.append(fake); fe.append(err)
                # Match the recipient mask to the exact usable mock support.
                oy.append(np.where(np.isfinite(err),f[t],np.nan)); oe.append(err)
                plans.append(plan)
                if target % 7 == 0:
                    for obj,y,se,model in ((real_quality,oy[-1],oe[-1],_model(oy[-1],oe[-1])),
                                           (mock_quality,fake,err,_model(fake,err))):
                        for k,val in zip(obj,_quality(y,se,model)):
                            obj[k].append(val)
            if not plans:
                continue
            totals["used_shots"]+=1
            indices=[p["row_index"] for p in plans]
            redshifts=cat.info["z_hetdex"][indices]
            for ys,es,stream,counts,dest in ((fy,fe,"nuisance_shuffle",mocked,leads),
                                            (oy,oe,"matched_observed",observed,real_leads)):
                summary,recs=screen_chunk(cat.wave,np.asarray(ys),np.asarray(es),redshifts,0,stream)
                counts.update(summary)
                for rec in recs:
                    design_row=plans[rec["row_index"]]
                    rec["row_index"]=design_row["row_index"]
                    rec["mock_design"]=design_row
                    dest.write(_json(rec)+"\n")
    denom=observed["spectra"]
    keys=("peaks","heii","ew","unknown")
    summary={"release":cat.release,"sampling":dict(totals),"expected_draw_fraction":1/STRATUM,
             "observed":dict(observed),"mock":dict(mocked),
             "rates_per_matched_spectrum":{
                 k:{"observed":observed[k]/denom if denom else None,"mock":mocked[k]/denom if denom else None}
                 for k in keys},"qa_subset":_quality_summary(real_quality,mock_quality),
             "block_channels":BLOCK,"minimum_donors":3,"maximum_donors":DONORS,
             "limitations":"Short blocks retain within-block covariance and sky-line wavelengths. Boundaries cut correlations. "
             "Genuine source lines can be transferred from donors. Rates condition on successful same-shot matching and shared masks."}
    (out/"nuisance-summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    return summary


def _detail_rows(cat,rows,y,e):
    from unknown_lines import unknown_spectrum
    for row in rows:
        i=row["row_index"]
        m=cat.metadata(i)
        diagnostics=spec.detail(cat.wave,y,e,row,density_unit=1e-17)
        if row["channel"]=="unknown":
            diagnostics["unknown_spectrum"]=unknown_spectrum(cat.wave,y,e,row,density_unit=1e-17)
        yield {**m,"catalog":m,**row,**diagnostics,"shuffled_detail_complete":True}


def records_for_vet(cat,out):
    from scan import jsonlines
    path=Path(out)/"nuisance-leads.jsonl.gz"
    for i,rows in groupby(jsonlines(path),lambda r:r["row_index"]):
        rows=list(rows)
        design=rows[0]["mock_design"]
        ids=[i,*design["donor_rows"]]
        f,e=spectra_for_rows(cat,ids)
        models=[_model(y,se) for y,se in zip(f,e)]
        y,err=_draw(models[0],models[1:],e[0],design["draw_seed"])
        yield from _detail_rows(cat,rows,y,err)


def raw_subsample_for_vet(cat,out,limit=1024):
    """Proportional stratified sample of raw-shuffled spectra, including zero leads."""
    from scan import BATCH,jsonlines
    strata=[]
    for start in range(0,cat.size,BATCH):
        n=min(BATCH,cat.size-start)
        j=start//BATCH
        if n>1 and j%8==random.Random(700001+j//8).randrange(8):
            strata.append((start,n))
    size=sum(n for _,n in strata)
    n_take=min(limit,size)
    phase=random.Random(SEED+97).random()
    sample=set()
    offset=0
    for start,n in strata:
        # Systematic rounding allocates proportional counts to blocks and
        # preserves the exact global budget. Uniformly sample within each.
        n_this=(int((offset+n)*n_take/size+phase)-int(offset*n_take/size+phase))
        sample.update(random.Random(SEED+97+start).sample(range(start,start+n),n_this))
        offset+=n
    rates=Counter()
    blocks=groupby(jsonlines(Path(out)/"shuffled.jsonl.gz"),lambda r:r['row_index']//BATCH*BATCH)
    for start,entries in blocks:
        subset=[r for r in entries if r["row_index"] in sample]
        if not subset:
            continue
        f,e=cat.read(start,min(start+BATCH,cat.size))
        y,se=spec.shuffle_fibres(f,e,seed=410087+start)
        for i,rows in groupby(subset,lambda r:r["row_index"]):
            rows=list(rows)
            rates.update(r["channel"] for r in rows)
            yield from _detail_rows(cat,rows,y[i-start],se[i-start])
    (Path(out)/"raw-shuffle-vet-sample.json").write_text(json.dumps({
        "seed":SEED+97,"raw_spectra":size,"selected_spectra":len(sample),
        "original_blocks":len(strata),"sampling":"proportional original raw-shuffle block strata",
        "selected_rows":sorted(sample),"sample_lead_counts":dict(rates)},indent=2)+"\n")
