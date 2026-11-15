"""Independent FEID event-epoch vs OMNI hourly Dst experiment.

FEID event onset_type encodes shock/SSC classification. Completed-event FD
amplitude is retrospective. Neither is a live neutron-monitor forecast.
Optional dependency: pyarrow (only required for Parquet input).
"""
from __future__ import annotations
import argparse
from bisect import bisect_left
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import statistics
from urllib.request import Request, urlopen

FEID = "https://huggingface.co/datasets/juliensimon/forbush-decreases/resolve/refs%2Fconvert%2Fparquet/default/train/0000.parquet"
OMNI = "https://huggingface.co/datasets/juliensimon/omni-solar-wind-parameters/resolve/refs%2Fconvert%2Fparquet/default/train/0000.parquet"
LIMIT_NT = -100.0
PREVIOUS_HOURS = 72
LAG_HOURS = 48
SEED = 20261115

def utc_hour(timestamp: datetime) -> int:
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)
    return int(timestamp.timestamp() // 3600)

def hour_to_month(hour: int) -> str:
    return datetime.fromtimestamp(hour * 3600, tz=timezone.utc).strftime("%Y-%m")

def fetch(url: str, path: Path) -> dict:
    req = Request(url, headers={"User-Agent": "Open-Space-research/0.2; public archive analysis"})
    path.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    with urlopen(req, timeout=180) as response, path.open("wb") as dest:
        while chunk := response.read(2 * 1024 * 1024):
            h.update(chunk)
            dest.write(chunk)
    return {"source": url, "sha256": h.hexdigest(), "bytes": path.stat().st_size}

def load_rows(path: Path, columns: list[str]) -> list[dict]:
    import pyarrow.parquet as pq  # optional, not needed for unit tests
    return pq.read_table(path, columns=columns).to_pylist()

def make_hourly_dst(rows: list[dict]) -> dict[int, float | None]:
    out = {}
    for row in rows:
        hour = utc_hour(row["datetime"])
        if hour in out:
            raise ValueError(f"OMNI has duplicate hourly timestamp: {hour}")
        value = row["dst_index_nt"]
        # NASA Dst fill=99999. Mirror usually maps this to null.
        valid = value is not None and -1500 < float(value) < 1500
        out[hour] = float(value) if valid else None
    if not out:
        raise ValueError("empty Dst input")
    return out

def severe_onsets_and_eligibility(dst: dict[int, float | None],
                                  threshold: float = LIMIT_NT,
                                  quiet_hours: int = PREVIOUS_HOURS) -> tuple[list[int], list[int]]:
    """Onset: first below-threshold hour after 72 *complete* quiet hours.

    Eligibility is pre-onset exposure with a valid, nonsevere current hour.
    Gaps or any previous storm hour veto both onset and eligibility.
    """
    window = deque()
    missing = storm = 0
    onsets, eligible = [], []
    for hour in range(min(dst), max(dst) + 1):
        value = dst.get(hour)
        if len(window) == quiet_hours and missing == 0 and storm == 0:
            if value is not None and value < threshold:
                onsets.append(hour)
            elif value is not None and value >= threshold:
                eligible.append(hour)
        if len(window) == quiet_hours:
            popped = window.popleft()
            missing -= (popped is None)
            storm -= (popped is not None and popped < threshold)
        window.append(value)
        missing += (value is None)
        storm += (value is not None and value < threshold)
    return onsets, eligible

def cluster_events(rows: list[dict], minimum_gap_hours: int = PREVIOUS_HOURS) -> list[dict]:
    """One earliest event per 72h episode, before reading outcome labels."""
    events = sorted((r for r in rows if r["datetime_utc"] is not None), key=lambda r: r["datetime_utc"])
    result = []
    for r in events:
        t = r["datetime_utc"]
        if not result or (t - result[-1]["last"]).total_seconds() >= minimum_gap_hours * 3600:
            result.append({"time": t, "last": t, "size": 1,
                           "onset_type": int(r["onset_type"]) if r["onset_type"] is not None else None})
        else:
            result[-1]["last"] = t
            result[-1]["size"] += 1
    return result

def lag_to_next(hour_float: float, onsets: list[int],
                horizon_hours: int = LAG_HOURS) -> float | None:
    i = bisect_left(onsets, hour_float)
    if i == len(onsets):
        return None
    lag = onsets[i] - hour_float
    return lag if 0 <= lag < horizon_hours else None

def wilson_bounds(hits: int, total: int, z: float = 1.96) -> list[float] | None:
    if total == 0:
        return None
    p = hits / total
    denom = 1 + z*z/total
    centre = (p + z*z/(2*total)) / denom
    half = z * ((p*(1-p)/total + z*z/(4*total*total)) ** 0.5) / denom
    return [round(max(0, centre-half), 6), round(min(1, centre+half), 6)]

def stratified_null(events: list[dict], onsets: list[int], hours_by_month: dict[str, list[int]],
                    repetitions: int = 2000, seed: int = SEED) -> dict:
    """Month/year matched random *eligible* UTC hours; episode-level draws."""
    rng = random.Random(seed)
    panels = [hours_by_month[hour_to_month(utc_hour(e["time"]))] for e in events]
    draws = [sum(lag_to_next(rng.choice(hours), onsets) is not None for hours in panels)
             for _ in range(repetitions)]
    draws.sort()
    return {"draws": repetitions, "seed": seed,
            "mean_hits": round(statistics.mean(draws), 5),
            "p025_hits": draws[int((repetitions-1)*0.025)],
            "p975_hits": draws[int((repetitions-1)*0.975)],
            "draws_sorted": draws}

def analyze(feid: list[dict], omni: list[dict], draws: int = 2000, seed: int = SEED) -> dict:
    dst = make_hourly_dst(omni)
    onsets, eligible = severe_onsets_and_eligibility(dst)
    eligible_set = set(eligible)
    by_month = defaultdict(list)
    for hour in eligible:
        by_month[hour_to_month(hour)].append(hour)
    episodes = cluster_events(feid)
    output = {
        "question": "FEID event epoch before independent OMNI Dst below -100 nT",
        "time_convention": "UTC; FEID exact minute; OMNI hourly timestamp",
        "dst_period": {"start": datetime.fromtimestamp(min(dst)*3600, timezone.utc).isoformat(),
                       "end": datetime.fromtimestamp(max(dst)*3600, timezone.utc).isoformat()},
        "dst_hours": len(dst), "dst_valid_hours": sum(v is not None for v in dst.values()),
        "severe_storm_onsets": len(onsets), "eligible_utc_hours": len(eligible),
        "raw_feid_rows": len(feid), "feid_72h_episodes": len(episodes),
        "filters": {"prior_quiet_hours": PREVIOUS_HOURS, "dst_threshold_nt_strict": LIMIT_NT,
                    "outcome_hours": LAG_HOURS, "first_event_per_72h_episode": True},
        "cohorts": {}
    }
    for name, (start, end) in {
        "development_1963_2000": (1963, 2000),
        "holdout_2001_2016": (2001, 2016),
        "all_1963_2016": (1963, 2016),
    }.items():
        subset = [e for e in episodes if start <= e["time"].year <= end]
        accepted = [e for e in subset if utc_hour(e["time"]) in eligible_set and
                    by_month.get(hour_to_month(utc_hour(e["time"])))]
        lags = [lag_to_next(e["time"].timestamp()/3600, onsets) for e in accepted]
        hit_lags = [x for x in lags if x is not None]
        bins = [sum(a <= t < b for t in hit_lags) for a,b in ((0,6),(6,12),(12,24),(24,48))]
        hits = len(hit_lags)
        null = stratified_null(accepted, onsets, by_month, draws,
                              seed + (1 if start==2001 else 0)) if accepted else None
        if null is not None:
            sorted_draws = null.pop("draws_sorted")
            null["greater_equal"] = sum(x >= hits for x in sorted_draws)
            null["empirical_upper_tail_p"] = round((1+null["greater_equal"])/(1+draws), 6)
            null["observed_to_null_mean_ratio"] = (
                round(hits/null["mean_hits"], 4) if null["mean_hits"] else None)
        type_counts = Counter(str(e["onset_type"]) for e in accepted)
        type_hits = Counter(str(e["onset_type"]) for e,t in zip(accepted,lags) if t is not None)
        # Label-stratified analyses make visible that the catalog codes contain
        # prior shock/SSC information. THESE ARE NOT REAL-TIME PREDICTORS.
        type_controls = {}
        for code in ("1", "9"):
            events_of_type = [e for e in accepted if str(e["onset_type"]) == code]
            if not events_of_type:
                continue
            observed = type_hits[code]
            controls = stratified_null(events_of_type, onsets, by_month, draws,
                                       seed + (100 if code == "1" else 900) +
                                       (1 if start == 2001 else 0))
            outcomes = controls.pop("draws_sorted")
            controls["hits_48h"] = observed
            controls["eligible_episodes"] = len(events_of_type)
            controls["upper_tail_p"] = round((1+sum(x >= observed for x in outcomes))/(1+draws), 6)
            controls["lower_tail_p"] = round((1+sum(x <= observed for x in outcomes))/(1+draws), 6)
            controls["observed_to_expected_ratio"] = (
                round(observed/controls["mean_hits"], 4) if controls["mean_hits"] else None)
            type_controls[code] = controls
        output["cohorts"][name] = {
            "episodes_in_period": len(subset), "eligible_episodes": len(accepted),
            "excluded_prior_storm_or_missing_dst": len(subset)-len(accepted),
            "hits_48h": hits, "rate": round(hits/len(accepted), 6) if accepted else None,
            "wilson_95_rate": wilson_bounds(hits, len(accepted)),
            "lag_bins_0_6_6_12_12_24_24_48": bins,
            "types_eligible": dict(type_counts), "types_hit": dict(type_hits),
            "type_specific_month_matched_controls": type_controls,
            "null_year_month_matched": null
        }
    return output

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feid", type=Path, default=Path("/tmp/forbush-feid.parquet"))
    parser.add_argument("--omni", type=Path, default=Path("/tmp/forbush-omni.parquet"))
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=2000)
    args = parser.parse_args()
    provenance = {}
    if args.fetch:
        provenance["feid"] = fetch(FEID, args.feid)
        provenance["omni"] = fetch(OMNI, args.omni)
    for name,path,uri in (("feid",args.feid,FEID),("omni",args.omni,OMNI)):
        if not path.is_file():
            parser.error(f"missing file {path}; pass --fetch or --{name} PATH")
        if name not in provenance:
            provenance[name] = {"source":uri,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                                "bytes":path.stat().st_size}
    feid = load_rows(args.feid, ["datetime_utc","onset_type"])
    omni = load_rows(args.omni, ["datetime","dst_index_nt"])
    result = analyze(feid,omni,draws=args.draws)
    result["source_files"] = provenance
    result["method"] = "examples/forbush_dst_independent.py; pyarrow for input only"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__ == "__main__":
    main()
