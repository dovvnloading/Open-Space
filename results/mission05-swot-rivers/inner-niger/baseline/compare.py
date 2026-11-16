"""Evaluate fixed calendar-window nadir-altimetry changes and archive denominators."""
import csv
import itertools
import json
import statistics
from datetime import date
from pathlib import Path

OUT=Path(__file__).resolve().parent


def read(target):
    with (OUT/target["file"]).open() as stream:
        return [(date.fromisoformat(r["time_utc"][:10]),float(r["water_level_m"]),float(r["uncertainty_m"]))
                for r in csv.DictReader(stream)]


def nearest(series,when,uncertainty):
    usable=[r for r in series if 0 < r[2] <= uncertainty and abs((r[0]-when).days)<=15]
    return min(usable,key=lambda r:(abs((r[0]-when).days),r[0])) if usable else None


def save(name,data):
    with (OUT/name).open("w",newline="") as stream:
        w=csv.DictWriter(stream,data[0],lineterminator="\n")
        w.writeheader()
        w.writerows(data)


def main():
    manifest=json.loads((OUT/"manifest.json").read_text())
    sites={}
    for site in ("mopti","dire"):
        sites[site]=sorted((t for t in manifest["targets"] if t["region"]==site
                            and t["first"]<"2019" and t["last"]>"2026-01-31"),
                           key=lambda t:t["distance_km"])
    residuals=[]
    results=[]
    for m,d in itertools.product(sites["mopti"],sites["dire"]):
        for unc in [0.30,0.15]:
            entries=[]
            series=[read(m),read(d)]
            for y in range(max(int(m["first"][:4]),int(d["first"][:4])),2026):
                row={"mopti_id":m["id"],"dire_id":d["id"],"primary_pair":m==sites["mopti"][0] and d==sites["dire"][0],
                     "uncertainty_limit_m":unc,"season":y}
                resolved=[]
                for label,s in zip(["mopti","dire"],series):
                    for month,when in [("july",date(y,7,22)),("january",date(y+1,1,26))]:
                        r=nearest(s,when,unc)
                        row[f"{label}_{month}_date"]=r[0].isoformat() if r else ""
                        row[f"{label}_{month}_m"]=r[1] if r else ""
                        row[f"{label}_{month}_u_m"]=r[2] if r else ""
                        resolved.append(r)
                complete=all(resolved)
                m_d=resolved[1][1]-resolved[0][1] if all(resolved[:2]) else None
                d_d=resolved[3][1]-resolved[2][1] if all(resolved[2:]) else None
                row["mopti_change_m"]=m_d if m_d is not None else ""
                row["dire_change_m"]=d_d if d_d is not None else ""
                row["dire_minus_mopti_change_m"]=d_d-m_d if complete else ""
                row["maximum_site_mismatch_days"]=max(abs((resolved[0][0]-resolved[2][0]).days),
                                                       abs((resolved[1][0]-resolved[3][0]).days)) if complete else ""
                entries.append(row)
            residuals.extend(entries)
            for stage_limit,site_skew in itertools.product([None,0.2,0.5],[None,7]):
                eligible=[r for r in entries if r["maximum_site_mismatch_days"]!=""
                          and (stage_limit is None or abs(r["mopti_change_m"])<=stage_limit)
                          and (site_skew is None or r["maximum_site_mismatch_days"]<=site_skew)]
                event=next((r for r in eligible if r["season"]==2024),None)
                pre=[r for r in eligible if r["season"]<=2021] # both visits predate Dec 2022 launch
                rank=1+sum(r["dire_change_m"]>=event["dire_change_m"] for r in pre) if event else ""
                record={"mopti_id":m["id"],"dire_id":d["id"],"primary_pair":entries[0]["primary_pair"],"u_limit_m":unc,
                        "mopti_stage_limit_m":stage_limit,"max_site_skew_d":site_skew,
                        "seasons_total":len(entries),"eligible_count":len(eligible),"eligible_preSWOT_count":len(pre),
                        "eligible_preSWOT_years":" ".join(str(r["season"]) for r in pre),
                        "target_season_eligible":event is not None,
                        "target_dire_change_m":event["dire_change_m"] if event else "",
                        "target_mopti_change_m":event["mopti_change_m"] if event else "",
                        "rank_vs_preSWOT_plus_event_desc":rank,
                        "comparison_count":len(pre)+1 if event else "",
                        "empirical_upper_tail_rank_fraction":rank/(len(pre)+1) if event else "",
                        "preSWOT_median_dire_change_m":statistics.median(r["dire_change_m"] for r in pre) if pre else "",
                        "preSWOT_min_dire_change_m":min((r["dire_change_m"] for r in pre),default=""),
                        "preSWOT_max_dire_change_m":max((r["dire_change_m"] for r in pre),default="")}
                results.append(record)
    save("seasonal_measurements.csv",residuals)
    save("seasonal_ranks.csv",results)


if __name__=="__main__":
    main()
