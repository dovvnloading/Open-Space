"""Count registered aperture and annulus photons in fixed energy/time bins of the discovery GTI."""
from collections import Counter
import csv
import gzip
import math
from pathlib import Path

import diagnostics
import ingest_primary as ingest
import numpy as np
from report import bounds

PHASES=("before","trigger","after")
REGIONS=("aperture","annulus")
BANDS=("soft","hard")
DELTAS=(180,600,1800,3600)
need=ingest.require


def intervals(lo,start,stop,hi):
    need(all(math.isfinite(t) for t in (lo,start,stop,hi)) and lo<=start<stop<=hi,"Invalid GTI or trial")
    a,b=start+(stop-start)/3,start+2*(stop-start)/3
    marks=[lo,*[start-d for d in reversed(DELTAS)],start,a,b,stop,*[stop+d for d in DELTAS],hi]
    clipped=[min(hi,max(lo,t)) for t in marks]
    names=("earlier_baseline","before_3600_1800","before_1800_600","before_600_180","before_180_0",
           "trial_first_third","trial_middle_third","trial_last_third","after_0_180","after_180_600",
           "after_600_1800","after_1800_3600","later_baseline")
    phases=["before"]*5+["trigger"]*3+["after"]*5
    return list(zip(names,phases,clipped[:-1],clipped[1:],strict=True))


def region_counts(counts,region):
    soft,hard=(counts[region,b] for b in BANDS)
    total=soft+hard
    lower,upper=bounds(hard,total)
    return {"soft":soft,"hard":hard,"total":total,"gross_hard_fraction":hard/total if total else None,
            "gross_hard_fraction_cp95":[lower,upper]}


def profile(rowfile,trigger):
    d=trigger["diagnostics"]
    need(ingest.sha(rowfile)==d["files"]["event_rows.csv.gz"],"Original event-row SHA-256 differs")
    lo,hi=map(float,d["gti_met_s"])
    start,stop=trigger["start_met_s"],trigger["stop_met_s"]
    ratio,duty=trigger["annulus_area_ratio"],trigger["live_fraction"]
    need(ratio>0 and 0<duty<=1 and math.isclose(start-lo,trigger["before_s"],abs_tol=1e-6,rel_tol=0)
         and math.isclose(hi-stop,trigger["after_s"],abs_tol=1e-6,rel_tol=0),"Original GTI, duty or annulus differs")
    tiles=intervals(lo,start,stop,hi)
    books=[Counter() for _ in tiles]
    totals={p:Counter() for p in PHASES}
    seen=set()
    peak=Counter()
    p0=trigger["peak_start_met_s"]
    p1=p0+trigger["peak_window_s"]
    with gzip.open(rowfile,"rt",encoding="ascii",newline="") as src:
        for r in csv.DictReader(src):
            idx=int(r["events_row_zero_based"])
            t,e=float(r["time_met_s"]),float(r["energy_eV"])
            region=r["region"]
            need(idx>=0 and idx not in seen and lo<=t<hi and 500<=e<=7000 and region in REGIONS and
                 int(r["grade"]) in {0,2,3,4,6},"Unexpected or repeated discovery event row")
            seen.add(idx)
            band="hard" if e>=2000 else "soft"
            phase="before" if t<start else "trigger" if t<stop else "after"
            need(phase==r["phase"] and int(r["in_peak"])==int(p0<=t<p1),"Published event phase differs")
            [k]=[i for i,(_,p,l,h) in enumerate(tiles) if p==phase and l<=t<h]
            books[k][region,band]+=1
            totals[phase][region,band]+=1
            if p0<=t<p1:
                peak[region]+=1
    checks={p:{r:region_counts(totals[p],r)["total"] for r in REGIONS} for p in PHASES}
    need(len(seen)==d["event_table_rows"] and all(checks[p]["aperture"]==trigger[k] for p,k in
         zip(PHASES,("before_counts","on_counts","after_counts"),strict=True)) and
         checks["trigger"]["annulus"]==trigger["annulus_counts"] and
         checks["before"]["annulus"]+checks["after"]["annulus"]==trigger["annulus_off_counts"] and
         peak["aperture"]==trigger["peak_counts"] and peak["annulus"]==trigger["peak_annulus_counts"],
         "Time/energy tally differs from original aperture, annulus or peak counts")
    def measure(counts):
        aperture=region_counts(counts,"aperture")
        annulus=region_counts(counts,"annulus")
        return {"aperture":aperture,"annulus":annulus,
                "scaled_background":{b:annulus[b]/ratio for b in BANDS},
                "net_aperture_counts":{b:aperture[b]-annulus[b]/ratio for b in BANDS}}
    return {"input_sha256":d["input_sha256"],"event_rows_sha256":ingest.sha(rowfile),"gti_met_s":[lo,hi],
            "trial_met_s":[start,stop],"bands_ev":{"soft":"500 <= energy < 2000","hard":"2000 <= energy <= 7000"},
            "source_center_sky_px":d["sky_center_px"],"radii_sky_px":d["apertures_sky_px"],"annulus_area_ratio":ratio,
            "live_fraction":duty,"source_event_rows":len(seen),"counts_by_phase":{p:measure(totals[p]) for p in PHASES},
            "bins":[{"name":name,"phase":p,"start_met_s":a,"stop_met_s":b,"wall_s":b-a,"live_s":(b-a)*duty,
                     **measure(c)} for (name,p,a,b),c in zip(tiles,books,strict=True)],
            "limits":"These are observed 0.5–7 keV aperture and annulus counts. The gross hard fraction conditions on "
            "the detected counts in each region; its two-sided 95% Clopper–Pearson interval uses the binomial counts. "
            "Intervals are per bin. Empty bins have no fraction or interval. The fractions include background and "
            "energy-dependent detector response. They cannot establish the source distance, absorption or a Galactic origin. "
            "The chosen field and binning are follow-up characterization; original instrument, 24-control and held-out results are unchanged."}


def regenerate(evt,event,expected,tmp):
    data,lo,hi=diagnostics._get_data(evt,event)
    need([lo,hi]==event["diagnostics"]["gti_met_s"],"Replayed GTI changed")
    dist=np.hypot(data["x"]-event["aperture_x"],data["y"]-event["aperture_y"])
    r=event["aperture_radius_px"]
    p0=event["peak_start_met_s"]
    path=Path(tmp)/"reproduced-event-rows.csv.gz"
    diagnostics._events_csv(path,data,dist<=r,(dist>=2*r)&(dist<=4*r),event["start_met_s"],event["stop_met_s"],
                            p0,p0+event["peak_window_s"])
    need(ingest.sha(path)==expected,"Pinned EVT2 cannot regenerate the original selected event rows")
    return expected


def markdown(data):
    lines=["# Time-resolved energy counts","",
           "Counts use the original source aperture and annulus. The hard band includes photons at 2.000 keV. "
           "The displayed hard fractions use the gross counts in each region.","",
           "| Bin | Live s | Source soft / hard | Source hard fraction (95%) | Annulus soft / hard | Annulus hard fraction (95%) |",
           "|---|---:|---:|---|---:|---|"]
    def fraction(r):
        if not r["total"]:
            return "—"
        low,high=r["gross_hard_fraction_cp95"]
        return f"{r['gross_hard_fraction']:.3f} ({low:.3f}–{high:.3f})"
    for b in data["bins"]:
        a,n=b["aperture"],b["annulus"]
        lines.append(f"| {b['name']} | {b['live_s']:.3f} | {a['soft']} / {a['hard']} | {fraction(a)} | "
                     f"{n['soft']} / {n['hard']} | {fraction(n)} |")
    lines+=["",data["limits"],"","[Counts, bin edges and phase sums](energy_bins.json) · [Original spatially selected events](event_rows.csv.gz)",""]
    return "\n".join(lines)
