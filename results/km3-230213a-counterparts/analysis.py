#!/usr/bin/env python3
"""One-event frozen KM3-230213A public-catalogue replication and bounded screen.

Standard library only. No failed request, unknown footprint, or unobserved
temporal interval is silently converted into a negative match.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import random
import re
import sys
import zipfile
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
EVENT = {"name": "KM3-230213A", "ra_deg": 94.3, "dec_deg": -7.8,
         "mjd_utc": 59988.0533299, "radii_deg": {"50": 1.2, "68": 1.5, "90": 2.2, "99": 3.0}}
T0 = datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(days=EVENT["mjd_utc"])
TNS_URL = "https://www.wis-tns.org/sites/default/files/tns_public_objects/tns_public_objects.csv.zip"
CHIME_URL = "https://assets.chime-frb.ca/catalogs/catalog2/releases/20251231/data/chimefrbcat2.csv"
N_SHIFT = 999
SEED = 230213
STATUS = "CONDITIONAL_CATALOGUE_ONLY"
receipts = {"run_utc": datetime.now(timezone.utc).isoformat(), "event": EVENT,
            "source_retrievals": [], "controls": {"seed": SEED, "valid_shifts": N_SHIFT,
            "method": "fixed event declination, deterministic random RA offsets 6.5..353.5 degrees, no reuse of target-selected known-17 list"}}

def sep(a, d, b, e):
    a, d, b, e = map(math.radians, (a, d, b, e))
    u = math.sin((d-e)/2)**2 + math.cos(d)*math.cos(e)*math.sin((a-b)/2)**2
    return 2*math.degrees(math.asin(min(1., math.sqrt(max(0.,u)))))

def parse_angle(value, is_ra=False):
    s = str(value or "").strip().replace("−", "-")
    if not s:
        return None
    if ":" not in s and " " not in s:
        try:
            x = float(s)
            return x if -90<=x<=360 else None
        except ValueError:
            return None
    tok = [t for t in re.split(r"[:\s]+", s) if t]
    if len(tok) < 2:
        return None
    try:
        a = list(map(float, tok[:3]))
    except ValueError:
        return None
    value = abs(a[0])+a[1]/60+(a[2]/3600 if len(a)>2 else 0)
    return math.copysign(value*(15 if is_ra else 1),a[0] if not s.startswith("-") else -1)

def parse_datetime(v):
    s = str(v or "").strip()
    if not s:
        return None
    try:
        t = datetime.fromisoformat(s.replace("Z","+00:00"))
        return t.replace(tzinfo=timezone.utc) if t.tzinfo is None else t.astimezone(timezone.utc)
    except ValueError:
        pass
    for pattern in ("%Y/%m/%d %H:%M:%S", "%Y-%m-%d", "%Y/%m/%d"):
        try:
            return datetime.strptime(s,pattern).replace(tzinfo=timezone.utc)
        except ValueError:
            pass
    return None

def pick(row,*keys):
    normalized = {re.sub("[^a-z0-9]","",str(k).lower()): v for k,v in row.items() if k}
    for key in keys:
        v=normalized.get(re.sub("[^a-z0-9]","",key.lower()))
        if v not in ("",None):
            return v
    return None

def download(label,url,timeout=115):
    record={"name":label,"url":url, "requested_utc": datetime.now(timezone.utc).isoformat()}
    try:
        req=Request(url,headers={"User-Agent":"Open-Space-KM3-public-research/1.0 (read-only, survey audit)","Accept":"*/*"})
        with urlopen(req,timeout=timeout) as res:
            # Cap to avoid failed mirrors returning enormous corrupt data.
            data=res.read(100_000_001)
            if len(data)>100_000_000:
                raise ValueError("exceeds 100 MB safety limit")
            record.update({"status":res.status,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),
                           "content_type":res.headers.get("Content-Type"),"final_url":res.url})
            return data
    except Exception as ex:
        record.update({"status":getattr(ex,"code",None),"error":str(ex)[:220],"outcome":"UNAVAILABLE_NOT_ZERO"})
        return None
    finally:
        receipts["source_retrievals"].append(record)

def known_replay():
    src=list(csv.DictReader((ROOT/"known_17.csv").open()))
    assert len(src)==17,"must transcribe all 17 published entries"
    out=[]
    for r in src:
        angle=sep(EVENT["ra_deg"],EVENT["dec_deg"],float(r["ra_deg"]),float(r["dec_deg"]))
        assert abs(angle-float(r["published_sep_deg"])) <=0.12,(r["name"],angle)
        assert angle<=3.0,(r["name"],angle)
        out.append({"name":r["name"],"known_id":r["id"],"distance_deg":round(angle,6),
           "published_rounded_distance_deg":r["published_sep_deg"],"p_chance":"NOT_ESTIMABLE_SELECTED_SAMPLE",
           "status":"KNOWN_ALREADY_PUBLISHED",**{"inside_"+k+"pct":angle<=v for k,v in EVENT["radii_deg"].items()}})
    write_csv("known_replay.csv",out,["known_id","name","distance_deg","published_rounded_distance_deg",
              "inside_50pct","inside_68pct","inside_90pct","inside_99pct","p_chance","status"])
    return {"count":len(out),"maximum_rounding_residual_deg":round(max(abs(r["distance_deg"]-float(r["published_rounded_distance_deg"])) for r in out),6),
            "inside_50pct":sum(r["inside_50pct"] for r in out),"inside_68pct":sum(r["inside_68pct"] for r in out),
            "inside_90pct":sum(r["inside_90pct"] for r in out),"inside_99pct":sum(r["inside_99pct"] for r in out)}

def write_csv(name,rows,cols):
    with (ROOT/name).open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=cols,extrasaction="ignore");w.writeheader();w.writerows(rows)

def tns_selection():
    data=download("TNS public objects snapshot",TNS_URL)
    if data is None:
        return None,{"status":"UNAVAILABLE","meaning":"cannot exclude an optical/TDE counterpart"}
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            paths=[x for x in z.namelist() if x.lower().endswith(".csv")]
            assert paths, "no CSV in TNS ZIP"
            with z.open(paths[0]) as f:
                lines=io.TextIOWrapper(f,encoding="utf-8-sig",errors="replace")
                # The TNS public ZIP starts with one timestamp preamble.
                head=next(lines,"")
                if "ra" in head.lower() and "," in head:
                    csvdata=io.StringIO(head+lines.read())
                else:
                    csvdata=io.StringIO(lines.read())
                reader=csv.DictReader(csvdata)
                rows=[]; total=0; parsed=0; eligible=0; fields=reader.fieldnames or []
                for r in reader:
                    total+=1
                    ra=parse_angle(pick(r,"ra","radeg"),True)
                    dec=parse_angle(pick(r,"declination","dec","decl","decdeg"))
                    when=parse_datetime(pick(r,"discoverydate","discoverydateutc","discovery_date","discovery"))
                    if ra is None or dec is None or when is None:
                        continue
                    parsed+=1
                    dt=(when-T0).total_seconds()/86400
                    typ=str(pick(r,"type","objecttype","classification") or "")
                    # Broad 30-d all types; 180-d TDE only, no posthoc expansion.
                    if abs(dt)>30 and not (abs(dt)<=180 and "TDE" in typ.upper()):
                        continue
                    eligible+=1
                    name=str(pick(r,"name","objname","objectname") or "")
                    prefix=str(pick(r,"nameprefix","prefix") or "AT")
                    if not name.startswith(("AT","SN")):
                        name=prefix+name
                    rows.append({"family":"TNS_OPTICAL" if abs(dt)<=30 else "TNS_TDE_LONG",
                                 "name":name,"ra":ra,"dec":dec,"dt_days":dt,
                                 "source_utc":when.isoformat(),"type":typ})
        if total < 50_000 or parsed < 40_000:
            raise ValueError("suspicious TNS table: rows=%d parsed=%d" % (total,parsed))
        return rows,{"status":"DOWNLOADED_FULL_PUBLIC_SNAPSHOT","raw_rows":total,
                     "rows_valid_position_and_discovery":parsed,"eligible_temporal_rows":eligible,
                     "columns":fields,"limitations":"TNS dates are discovery dates, not onset/peak; public classification and survey coverage incomplete"}
    except Exception as ex:
        return None,{"status":"PARSE_UNAVAILABLE","error":str(ex)[:250],"meaning":"cannot exclude a TNS transient"}

def chime_selection():
    data=download("CHIME/FRB Catalog 2 20251231",CHIME_URL)
    if data is None:
        return None,{"status":"UNAVAILABLE","meaning":"cannot exclude FRB candidate"}
    try:
        rows=[]; seen=set()
        reader=csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
        fields=reader.fieldnames or []
        total=0; parsed=0
        min_mjd=float("inf");max_mjd=-float("inf");nearest_dt=float("inf")
        near_7d=0;near_30d=0
        for r in reader:
            total+=1
            name=str(pick(r,"tns_name","event_id","tnsname") or "")
            ra=parse_angle(pick(r,"ra","ra_deg"),False)
            dec=parse_angle(pick(r,"dec","dec_deg"),False)
            mjd=pick(r,"mjd_inf","mjd_400","mjd")
            try: t=float(mjd)
            except (ValueError,TypeError):continue
            if ra is None or dec is None or not (40000<t<100000):continue
            parsed+=1
            min_mjd=min(min_mjd,t);max_mjd=max(max_mjd,t)
            nearest_dt=min(nearest_dt,abs(t-EVENT["mjd_utc"]))
            near_7d+=abs(t-EVENT["mjd_utc"])<=7
            near_30d+=abs(t-EVENT["mjd_utc"])<=30
            if not name: name="CHIME-%d" % total
            if name in seen: continue
            seen.add(name)
            delta=t-EVENT["mjd_utc"]
            if abs(delta)>1:continue
            rows.append({"family":"CHIME_FRB","name":name,"ra":ra,"dec":dec,
                "dt_days":delta,"source_utc":"MJD_UTC_"+str(mjd),"type":"FRB",
                "ra_err":pick(r,"ra_err","ra_err_deg"),"dec_err":pick(r,"dec_err","dec_err_deg")})
        if total < 3000 or parsed < 2500:
            raise ValueError("suspicious CHIME Catalog 2: rows=%d parsed=%d" % (total,parsed))
        return rows,{"status":"DOWNLOADED_CATALOG_2","raw_rows":total,"parsed_rows":parsed,
                    "min_mjd":min_mjd,"max_mjd":max_mjd,
                    "nearest_event_abs_days":round(nearest_dt,7),
                    "subbursts_within_7days":near_7d,
                    "subbursts_within_30days":near_30d,
                    "eligible_temporal_unique_events":len(rows),"columns":fields,
                    "limitations":"CHIME exposure nonuniform; centroid-only 3-degree membership is not localization overlap and multi-beam errors require native region"}
    except Exception as ex:
        return None,{"status":"PARSE_UNAVAILABLE","error":str(ex)[:250],"meaning":"cannot exclude an FRB"}

def analyze_family(family,rows,shifts):
    # Returns conditional empirical sky-shift calibration, not astrophysical posterior.
    area=2*math.pi*(1-math.cos(math.radians(3)))*(180/math.pi)**2
    seen=[]
    for r in rows:
        d=sep(EVENT["ra_deg"],EVENT["dec_deg"],r["ra"],r["dec"])
        if d<=3:
            seen.append((r,d))
    shift_counts=[]
    for shift in shifts:
        ra=(EVENT["ra_deg"]+shift)%360
        shift_counts.append(sum(sep(ra,EVENT["dec_deg"],r["ra"],r["dec"])<=3 for r in rows))
    obs=len(seen)
    pcount=(1+sum(c>=obs for c in shift_counts))/(len(shift_counts)+1)
    plower=(1+sum(c<=obs for c in shift_counts))/(len(shift_counts)+1)
    # For each row: "as good or better" spatial AND temporal match, within pre-fixed window,
    # and *any* object of the full same temporal-selected catalogue.
    results=[]
    for r,d in seen:
        rivals=0
        for shift in shifts:
            ra=(EVENT["ra_deg"]+shift)%360
            if any(sep(ra,EVENT["dec_deg"],s["ra"],s["dec"])<=d and
                   abs(s["dt_days"])<=abs(r["dt_days"])+1e-12 for s in rows):
                rivals+=1
        p=(1+rivals)/(len(shifts)+1)
        results.append({"family":family,"name":r["name"],"event_dt_days":round(r["dt_days"],8),
         "sep_deg":round(d,7),"within_1h":abs(r["dt_days"])<=1/24,
         "within_90pct":d<=2.2,"within_99pct":d<=3.0,
         "ra_deg":r["ra"],"dec_deg":r["dec"],"date_from_catalog":r["source_utc"],
         "source_type":r.get("type",""),"status":"UNVERIFIED_CATALOGUE_COINCIDENCE",
         "p_conditional_skyshift":round(p,6),
         "p_upper_bound_16_trials":round(min(1.,16*p),6),
         "p_physical_association":"NOT_ESTIMABLE_NO_EXPOSURE_FUNCTION",
         "localization_warning":"centroid-only" if family=="CHIME_FRB" else "discovery date, not onset",
         "ra_err":r.get("ra_err",""),"dec_err":r.get("dec_err","")})
    mean=sum(shift_counts)/len(shift_counts)
    return results,{"status":"CONDITIONAL_SCREENER_ONLY", "temporal_parent_count":len(rows),
      "observed_in_3deg":obs, "empirical_mean_count":round(mean,5),
      "p_shift_count_ge_observed":round(pcount,6),
      "p_shift_count_le_observed":round(plower,6),
      "skyshift_min_count":min(shift_counts),"skyshift_max_count":max(shift_counts),
      "skyshift_n_zero":sum(n==0 for n in shift_counts),
      "uniform_sky_density_per_deg2":round(len(rows)/41252.96125,8),
      "uniform_sky_expected_3deg":round(len(rows)*area/41252.96125,6),
      "uniform_density_caveat":"diagnostic only, not source or chance association posterior",
      "n_shifts":len(shifts),"limitations":"RA-shift ignores survey pointing, cadence and Galactic gradients; no complete selection function"}

def main():
    base=known_replay()
    rr=random.Random(SEED)
    shifts=[rr.uniform(6.5,353.5) for i in range(N_SHIFT)]
    # This fixes the shift set before reading either new dataset.
    tns, ts=tns_selection()
    chime,cs=chime_selection()
    summaries={"event":EVENT, "known_reproduction":base,
       "public_catalogue_status":{"TNS":ts,"CHIME":cs},
       "statistical_interpretation":"conditional sky-shift p is a catalog-only diagnostic, not source association posterior; no globally validated p for incomplete search families",
       "families":{}}
    new=[]
    if tns is not None:
        for family in ("TNS_OPTICAL","TNS_TDE_LONG"):
            matches,stats=analyze_family(family,[r for r in tns if r["family"]==family],shifts)
            summaries["families"][family]=stats;new.extend(matches)
    if chime is not None:
        matches,stats=analyze_family("CHIME_FRB",chime,shifts)
        summaries["families"]["CHIME_FRB"]=stats;new.extend(matches)
    # Freeze all 2023-epoch temporally eligible parent objects so another
    # implementation can independently verify the positions and sky-shift null.
    parent=(tns or [])+(chime or [])
    write_csv("temporal_parent_sample.csv",parent,
       ["family","name","ra","dec","dt_days","source_utc","type","ra_err","dec_err"])
    new.sort(key=lambda r:(r["p_conditional_skyshift"],r["event_dt_days"],r["name"]))
    write_csv("new_associations.csv",new,["family","name","event_dt_days","sep_deg",
      "within_1h","within_90pct","within_99pct","ra_deg","dec_deg","date_from_catalog",
      "source_type","status","p_conditional_skyshift","p_upper_bound_16_trials",
      "p_physical_association","localization_warning","ra_err","dec_err"])
    (ROOT/"control_statistics.json").write_text(json.dumps(summaries,indent=2,sort_keys=True)+"\n")
    (ROOT/"retrieval_receipts.json").write_text(json.dumps(receipts,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"known":base, "families":{k:{"status":v["status"],"observed":v["observed_in_3deg"],
     "p_shift_count":v["p_shift_count_ge_observed"]} for k,v in summaries["families"].items()},
     "download_status":{k:v["status"] for k,v in summaries["public_catalogue_status"].items()},
     "new_catalog_coincidences":len(new)},sort_keys=True))
    # Failed archives are a real audit outcome, NOT a negative match.
    return 0

if __name__=="__main__":
    sys.exit(main())
