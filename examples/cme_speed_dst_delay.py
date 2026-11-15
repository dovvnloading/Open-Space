"""Cross-cycle CME speed vs shock-to-Dst-minimum latency archival test.

Scope is *conditioned on storms already classified as severe*. The source
catalog's CME association is retrospective. CME sky-plane speed can be known
pre-arrival, but this retrospective association is NOT an operational forecast.
The 700/1000 km/s cutpoints are exploratory; paper novelty is NOT asserted.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import random
import statistics
from urllib.request import Request,urlopen

SOURCE_COMMIT="97747daa2ab3256986bdeafa707dd8214da71911"
URL="https://raw.githubusercontent.com/vema1818/GeoMagStorms_SolarSources/"+SOURCE_COMMIT+"/List_GMS_Storms.csv"
SEED=20261115

def parse_utc(s: str) -> datetime | None:
    try:
        return datetime.strptime(s.strip(),"%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    except (ValueError,AttributeError):
        return None

def parse_catalog(raw: str) -> tuple[list[dict],dict]:
    records=list(csv.DictReader(io.StringIO(raw)))
    events=[]
    counts=Counter()
    for r in records:
        peak=parse_utc(r["Peak Date/Time"])
        shock=parse_utc(r["icme onsett (yyyy-mm-dd hh:mm)"])
        if peak is None:
            counts["invalid_peak"]+=1
            continue
        if shock is None:
            counts["no_shock_timestamp"]+=1
            continue
        h=(peak-shock).total_seconds()/3600
        if not 0<=h<100:
            counts["invalid_arrival_to_minimum_delay"]+=1
            continue
        counts["icme_with_shock"]+=1
        try:
            cme=float(r["Cme speed (km/s)"])
        except (ValueError,TypeError):
            counts["missing_cme_speed"]+=1
            continue
        if not 0<cme<6000:
            counts["invalid_cme_speed"]+=1
            continue
        events.append({
            "event_id":int(r["No"]),"year":peak.year,
            "shock_epoch":shock.timestamp(),"delay_hours":h,
            "cme_speed_km_s":cme,"dst_min_nt":float(r["DST Peak (nT)"]),
            "bz_min_nt":float(r["Bz Peak(nT)"]),
            "quiet_region":int("quiet" in r["ACTIVE REGION"].lower())
        })
    counts["total_storm_rows"]=len(records)
    counts["cme_speed_analyzable"]=len(events)
    return events,dict(counts)

def pearson(pairs: list[tuple[float,float]]) -> float:
    n=len(pairs)
    if n<3: raise ValueError("need at least 3 pairs")
    sx=sum(x for x,y in pairs)
    sy=sum(y for x,y in pairs)
    sxx=sum(x*x for x,y in pairs)
    syy=sum(y*y for x,y in pairs)
    sxy=sum(x*y for x,y in pairs)
    den=(sxx-sx*sx/n)*(syy-sy*sy/n)
    if den<=0: raise ValueError("constant input")
    return (sxy-sx*sy/n)/math.sqrt(den)

def correlation(a: list[dict]) -> float:
    return pearson([(v["cme_speed_km_s"],v["delay_hours"]) for v in a])

def median_gap(a: list[dict]) -> float | None:
    slow=[x["delay_hours"] for x in a if x["cme_speed_km_s"]<700]
    fast=[x["delay_hours"] for x in a if x["cme_speed_km_s"]>=1000]
    if not slow or not fast:return None
    return statistics.median(slow)-statistics.median(fast)

def cohort(a: list[dict]) -> dict:
    slow=[x["delay_hours"] for x in a if x["cme_speed_km_s"]<700]
    fast=[x["delay_hours"] for x in a if x["cme_speed_km_s"]>=1000]
    return {"n":len(a),"pearson_r":round(correlation(a),6),
            "slow_cme_lt700_n":len(slow),"slow_median_h":round(statistics.median(slow),6) if slow else None,
            "fast_cme_ge1000_n":len(fast),"fast_median_h":round(statistics.median(fast),6) if fast else None,
            "slow_minus_fast_median_h":round(median_gap(a),6) if slow and fast else None}

def quantiles(a: list[float]) -> list[float]:
    a=sorted(a)
    return [round(a[int((len(a)-1)*v)],6) for v in (.025,.975)]

def resamples(a: list[dict], draws: int, seed: int) -> dict:
    rand=random.Random(seed)
    r_samples,gap_samples=[],[]
    for _ in range(draws):
        sample=rand.choices(a,k=len(a))
        r_samples.append(correlation(sample))
        gap=median_gap(sample)
        if gap is not None:gap_samples.append(gap)
    years=sorted(set(x["year"] for x in a))
    per_year={yr:[x for x in a if x["year"]==yr] for yr in years}
    year_r=[]
    for _ in range(draws):
        sampled_years=rand.choices(years,k=len(years))
        sample=[item for yr in sampled_years for item in per_year[yr]]
        year_r.append(correlation(sample))
    return {"draws":draws,"seed":seed,
            "pearson_95_percentile_interval":quantiles(r_samples),
            "median_gap_hours_95_percentile_interval":quantiles(gap_samples),
            "calendar_year_block_pearson_95_interval":quantiles(year_r)}

def era(year: int) -> int:
    return 1996+5*((year-1996)//5)

def era_permutation(a: list[dict], permutations: int, seed: int) -> dict:
    """Shuffle delays within five-year blocks, not across solar cycles."""
    rand=random.Random(seed)
    original=correlation(a)
    groups=defaultdict(list)
    for i,x in enumerate(a):groups[era(x["year"])].append(i)
    pairs=[(x["cme_speed_km_s"],x["delay_hours"]) for x in a]
    null=[]
    for _ in range(permutations):
        copied=[list(p) for p in pairs]
        for group in groups.values():
            vals=[copied[i][1] for i in group]
            rand.shuffle(vals)
            for i,new in zip(group,vals):copied[i][1]=new
        null.append(pearson(copied))
    return {"permutations":permutations,"seed":seed,
            "null_mean_r":round(statistics.mean(null),6),"null_95_r":quantiles(null),
            "two_sided_empirical_p":round((1+sum(abs(x)>=abs(original) for x in null))/(permutations+1),7)}

def solve(a: list[list[float]], b: list[float]) -> list[float]:
    n=len(b)
    augmented=[list(a[i])+[b[i]] for i in range(n)]
    for k in range(n):
        j=max(range(k,n),key=lambda row:abs(augmented[row][k]))
        augmented[k],augmented[j]=augmented[j],augmented[k]
        pivot=augmented[k][k]
        if abs(pivot)<1e-9:raise ValueError("singular design")
        augmented[k]=[x/pivot for x in augmented[k]]
        for i in range(n):
            if i!=k:
                scale=augmented[i][k]
                augmented[i]=[v-scale*w for v,w in zip(augmented[i],augmented[k])]
    return [row[-1] for row in augmented]

def adjusted_partial_correlation(a: list[dict]) -> float:
    """Diagnostic-only adjustment for Dst/Bz observed AFTER shock and epoch."""
    periods=sorted(set(era(v["year"]) for v in a))[1:]
    X=[[1.,v["dst_min_nt"]/100.,v["bz_min_nt"]/10.,float(v["quiet_region"])] +
       [float(era(v["year"])==p) for p in periods] for v in a]
    k=len(X[0])
    XX=[[sum(x[i]*x[j] for x in X) for j in range(k)] for i in range(k)]
    def residual(target):
        vals=[v[target] for v in a]
        Xy=[sum(x[i]*y for x,y in zip(X,vals)) for i in range(k)]
        beta=solve(XX,Xy)
        return [y-sum(b*z for b,z in zip(beta,x)) for x,y in zip(X,vals)]
    return pearson(list(zip(residual("cme_speed_km_s"),residual("delay_hours"))))

def independent_shock_episodes(a: list[dict], min_gap_hours: float = 72) -> list[dict]:
    """Retain first shock after the prior retained shock, for timing sensitivity."""
    out=[]
    for item in sorted(a,key=lambda x:x["shock_epoch"]):
        if not out or item["shock_epoch"]-out[-1]["shock_epoch"]>=min_gap_hours*3600:
            out.append(item)
    return out

def analyze(a: list[dict], counts: dict, bootstrap: int = 4000,
            permutations: int = 10000) -> dict:
    if len(a)<3:raise ValueError("not enough annotated storms")
    pre=[x for x in a if x["year"]<=2010]
    post=[x for x in a if x["year"]>=2011]
    iso=independent_shock_episodes(a)
    trimmed=[x for x in a if 100<=x["cme_speed_km_s"]<=2500 and x["delay_hours"]<49]
    return {
        "claim_class":"exploratory_archive_association_not_prospective_prediction",
        "sample_definition":"Severe Dst<=-100 nT catalog, 1996-2025, ICME shock at L1 and valid CME plane-of-sky speed",
        "source_url":URL,"source_commit":SOURCE_COMMIT,
        "time_zone":"UTC","cutoffs_km_s":{"slow_strictly_below":700,"fast_at_least":1000},
        "source_counts":counts,
        "cohorts":{"all":cohort(a),"earlier_1996_2010":cohort(pre),
                   "later_2011_2025":cohort(post),
                   "first_shock_per_72h":cohort(iso),"plausible_speed_delay_trim":cohort(trimmed)},
        "uncertainty":resamples(a,bootstrap,SEED),
        "five_year_stratified_permutation":era_permutation(a,permutations,SEED+1),
        "partial_r_dst_bz_quiet_source_era_diagnostic_only":round(adjusted_partial_correlation(a),6),
        "limitations":["The source catalog selects on Dst severity; it cannot predict WHETHER a severe storm occurs.",
                       "CME/storm matching is made retrospectively by observers.",
                       "1996-2010 vs 2011-2025 are retrospective era comparisons, not sealed blind holdouts.",
                       "Speed cutpoints were selected after inspecting this dataset; no confirmatory threshold test.",
                       "Dst and Bz are end-of-event diagnostic adjustment covariates, NOT usable pre-arrival.",
                       "Published literature already establishes storm timing and shock/sheath physics; no discovery priority asserted."]
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch",action="store_true")
    p.add_argument("--input",type=Path,default=Path("/tmp/geo-storms.csv"))
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=4000)
    p.add_argument("--permutations",type=int,default=10000)
    args=p.parse_args()
    if args.fetch:
        req=Request(URL,headers={"User-Agent":"Open-Space-research/0.2"})
        with urlopen(req,timeout=90) as src,args.input.open("wb") as dest:
            while chunk:=src.read(1024*1024):dest.write(chunk)
    if not args.input.is_file():p.error("missing input; pass --fetch or --input")
    blob=args.input.read_bytes()
    events,counts=parse_catalog(blob.decode("utf-8-sig"))
    result=analyze(events,counts,bootstrap=args.bootstrap,permutations=args.permutations)
    result["source_sha256"]=hashlib.sha256(blob).hexdigest()
    result["source_bytes"]=len(blob)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
