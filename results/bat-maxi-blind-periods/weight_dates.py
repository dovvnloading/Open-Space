"""Audit the calendar dates of extreme formal weights in frozen daily products."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from field_checks import load_archives, read_rows
from search import bounded_uncertainties, kish

MIN_ROWS, MIN_SPAN = 40, 60.
CONTROL_IDS = {"VelaX-1","CenX-3","HerX-1"}
CHECK_DAYS = (58571,58572)


def summarize(values):
    return {"min":float(np.min(values)), "median":float(np.median(values)),
            "max":float(np.max(values))} if len(values) else None


def number(values, fn=np.median):
    return float(fn(values)) if len(values) else None


def write_table(path, rows):
    if not rows:
        return
    opener = gzip.open if path.suffix==".gz" else open
    with opener(path,"wt",newline="") as out:
        writer = csv.DictWriter(out,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def band_audit(root, output, band, roster, receipts):
    wanted = [sid for sid in roster if band=="bat" or roster[sid].get("maxi_id")]
    valid = [receipts[sid][band] for sid in wanted if receipts.get(sid,{}).get(band,{}).get("status")=="ok"]
    if not valid:
        raise ValueError("No usable daily archive receipts for "+band)
    first = int(np.floor(min(r["first_mjd"] for r in valid)))
    last = int(np.floor(max(r["last_mjd"] for r in valid)))
    days = np.arange(first,last+1)
    shape = (len(days),len(wanted))
    # A missing value denotes absence of an observation on that source-day.
    share = np.full(shape,np.nan,dtype=np.float64)
    capped_share = np.full_like(share,np.nan)
    error_ratio = np.full_like(share,np.nan)
    tied = np.zeros(len(days))
    max_holders = [[] for _ in days]
    pool_logs, rows, controls = [], [], []
    rejected = Counter()
    seen = set()
    by_shard = {}
    for source in wanted:
        part = receipts[source]["part"]
        by_shard.setdefault(part,[]).append(source)
    for part,ids in sorted(by_shard.items()):
        archive = root/f"curves-{part:02}.tar.gz"
        curves = load_archives([archive],{f"{sid}.{band}" for sid in ids})
        for sid in ids:
            product = receipts[sid][band]
            curve = curves.get(f"{sid}.{band}")
            if curve is None:
                if product.get("status")=="ok":
                    raise ValueError(f"Missing frozen {band} curve for successful receipt {sid}")
                rejected["product unavailable"]+=1
                continue
            t,_,e = curve
            good = np.isfinite(t)&np.isfinite(e)&(e>0)
            t,e = t[good],e[good]
            order = np.argsort(t)
            t,e = t[order],e[order]
            if len(t)<MIN_ROWS or np.ptp(t)<MIN_SPAN:
                rejected["fewer than 40 epochs or less than 60 days"]+=1
                continue
            if sid in seen:
                raise ValueError("Duplicate source in frozen audit: "+sid)
            seen.add(sid)
            capped,_ = bounded_uncertainties(t,e)
            raw = (e.min()/e)**2
            new = (capped.min()/capped)**2
            raw_total = raw.sum()
            # Keep the total native inverse variance in log units until the
            # final pooling across different sources within an instrument.
            pool_logs.append(float(np.log(raw_total)-2*np.log(e.min())))
            raw,new = raw/raw_total,new/new.sum()
            dates = np.floor(t).astype(int)
            d,starts = np.unique(dates,return_index=True)
            i = len(rows)
            share[d-first,i] = np.add.reduceat(raw,starts)
            capped_share[d-first,i] = np.add.reduceat(new,starts)
            error_ratio[d-first,i] = np.minimum.reduceat(e,starts)/np.median(e)
            maxima = np.unique(dates[e==e.min()])
            for day in maxima:
                tied[day-first] += 1/len(maxima)
                max_holders[day-first].append(i)
            per = {"instrument":band,"source_id":sid,"name":roster[sid]["name"],
                   "sha256_archive_product":product["sha256"], "samples":len(e),
                   "span_days":float(np.ptp(t)),"median_quoted_error":float(np.median(e)),
                   "minimum_quoted_error":float(e.min()),"minimum_over_median_error":float(e.min()/np.median(e)),
                   "raw_kish":kish(e),"bounded_kish":kish(capped),"tied_max_weight_dates":len(maxima),
                   "first_max_weight_mjd":int(maxima[0]),"last_max_weight_mjd":int(maxima[-1]),
                   "raw_max_observation_share":float(raw.max()),
                   "raw_max_day_share":float(share[d-first,i].max()),
                   "bounded_max_day_share":float(capped_share[d-first,i].max())}
            rows.append(per)
            if sid in CONTROL_IDS:
                controls.append({**per, "check_dates":[{"mjd":v,"observed":bool(v in d),
                    "raw_weight_share":float(share[v-first,i]) if v in d else None,
                    "bounded_weight_share":float(capped_share[v-first,i]) if v in d else None,
                    "quoted_error_over_source_median":float(error_ratio[v-first,i]) if v in d else None}
                    for v in CHECK_DAYS]})
    n = len(rows)
    if not n:
        raise ValueError("No eligible source records for "+band)
    pool = np.exp(np.array(pool_logs)-max(pool_logs))
    pool /= pool.sum()
    date_rows = []
    for idx,day in enumerate(days):
        raw = share[idx,:n]
        observed = np.isfinite(raw)
        if not observed.any():
            continue
        nw = int(observed.sum())
        owners = np.asarray(max_holders[idx],dtype=int)
        peak_raw = raw[owners]
        bounded = capped_share[idx,:n]
        ratios = error_ratio[idx,:n]
        date_rows.append({"instrument":band,"mjd":int(day),"eligible_sources":n,"observed_sources":nw,
            "sources_with_max_weight":len(owners),"fraction_eligible_with_max":len(owners)/n,
            "fraction_observed_with_max":len(owners)/nw,"tie_adjusted_source_count":float(tied[idx]),
            "tie_adjusted_fraction_eligible":float(tied[idx]/n),
            "source_normalized_weight_fraction_all":float(np.nansum(raw,dtype=float)/n),
            "source_normalized_bounded_fraction_all":float(np.nansum(bounded,dtype=float)/n),
            "native_pooled_weight_fraction_all":float(np.sum(raw[observed]*pool[observed])),
            "median_weight_fraction_observed":float(np.median(raw[observed])),
            "max_weight_fraction_observed":float(np.max(raw[observed])),
            "median_error_over_source_median_observed":float(np.median(ratios[observed])),
            "median_weight_fraction_max_holders":number(peak_raw),
            "median_bounded_weight_fraction_max_holders":number(bounded[owners]),
            "median_error_over_source_median_max_holders":number(ratios[owners]),
            "max_holders_with_25pct_or_more":int(np.sum(peak_raw>=.25)),
            "observed_with_25pct_or_more":int(np.sum(raw[observed]>=.25)),
            "observed_with_50pct_or_more":int(np.sum(raw[observed]>=.5))})
    rank = sorted(date_rows,key=lambda r:(-r["tie_adjusted_source_count"],r["mjd"]))
    write_table(output/f"{band}-weight-dates.csv.gz",date_rows)
    write_table(output/f"{band}-source-weights.csv.gz",rows)
    selected = [day-first for day in CHECK_DAYS if first<=day<=last]
    union = np.unique([i for d in selected for i in max_holders[d]]).astype(int)
    observed_either = np.any(np.isfinite(share[selected,:n]),axis=0)
    total_group = np.nansum(share[selected,:n],axis=0)
    total_bounded = np.nansum(capped_share[selected,:n],axis=0)
    cluster = {"dates":CHECK_DAYS,"eligible_sources":n,"observed_on_either":int(observed_either.sum()),
               "unique_sources_with_max_on_either":len(union),"fraction_eligible":len(union)/n,
               "raw_share_for_max_holders":summarize(total_group[union]),
               "bounded_share_for_max_holders":summarize(total_bounded[union]),
               "source_normalized_raw_weight_fraction_all":float(total_group.sum()/n),
               "native_pooled_weight_fraction_all":float(total_group@pool)}
    return {"instrument":band,"roster_products":len(wanted),"eligible_sources":n,"excluded":rejected,
            "sources_with_ties_across_dates":sum(r["tied_max_weight_dates"]>1 for r in rows),
            "source_maximum_daily_fraction":summarize([r["raw_max_day_share"] for r in rows]),
            "source_maximum_daily_fraction_bounded":summarize([r["bounded_max_day_share"] for r in rows]),
            "top_dates_by_tie_adjusted_source_count":rank[:10],
            "specified_dates":[r for r in date_rows if r["mjd"] in CHECK_DAYS], "specified_two_day_union":cluster,
            "named_examples":controls}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    roster={r["id"]:r for r in json.loads((args.evidence/"inventory.json").read_text())["sources"]}
    receipts={}
    for path in sorted(args.evidence.glob("ingest-*.jsonl")):
        part=int(path.stem.split("-")[-1])
        for row in read_rows(path):
            if row["id"] in receipts:
                raise ValueError("Duplicate archive receipt for "+row["id"])
            receipts[row["id"]]={**row,"part":part}
    if set(receipts)!=set(roster):
        raise ValueError("Archive receipt census differs from the BAT roster")
    metadata={"eligibility":"At least 40 actual valid epochs and at least 60 elapsed days. Original inverse-variance concentration and source flux do not enter eligibility. MAXI comparison uses only uniquely associated BAT-roster sources.",
        "dates":"Each MJD day is the half-open interval from its integer value through the next integer. At tied minima each affected date is counted and receives 1/k votes in the tie-adjusted ranking. Per-day weight shares sum all retained observations within that day.",
        "weight":"Raw weights are inverse quoted variance. For each source all weights sum to one. The mean over every eligible source includes zero when it was unobserved; the separate native pooled fraction instead sums inverse variances across the instrument's sources before normalization. Medians labelled observed exclude absent sources. The error ratio is the lowest quoted error on the source-day divided by its source's median quoted error.",
        "controls":"The second-pass error cap provides a deterministic leverage comparison. The matched MAXI curves provide a different cadence and uncertainty pipeline. Vela X-1, Cen X-3 and Her X-1 are named examples selected after the original investigation. Report the complete date ranking to distinguish this selection from the all-source census.",
        "limits":"Multiple sources can share exposures, pointings, detector conditions and processing. These empirical concentrations do not identify the cause or estimate a false-alarm probability for an archive defect; they do not establish an error in any quoted uncertainty or an astrophysical event.",
        "bands":[band_audit(args.evidence,args.output,b,roster,receipts) for b in ("bat","maxi")]}
    (args.output/"weight-date-summary.json").write_text(json.dumps(metadata,indent=2,allow_nan=False)+"\n")
    for b in metadata["bands"]:
        top=b["top_dates_by_tie_adjusted_source_count"][0]
        print(b["instrument"], b["eligible_sources"], "eligible; most common maximum",top["mjd"],
              top["sources_with_max_weight"],"sources",flush=True)


if __name__=="__main__":
    main()
