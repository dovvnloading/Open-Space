"""Cross-match all catalogue local-time predictions with Fermi-GBM event UTCs."""
import argparse
import csv
from collections import Counter
from datetime import datetime
from pathlib import Path

from lxml import html
import numpy as np
from astropy.time import Time, TimeDelta

HERE = Path(__file__).resolve().parent
SHIFTS = (-69.184, 0.0, 69.184, 138.368)
ALIASES = {"Swift J1555.2-5402":"SGR J1555.2-5402",
           "Swift J1818.0-1607":"SGR J1818.0-1607",
           "PSR J1846-0258":"SGR J1846-0258","4U 0142+61":"AXP 4U 0142+61",
           "Swift J1822.3-1606":"SGR J1822.3-1606"}


def read_gbm(path):
    events=[]
    for tr in html.parse(str(path)).xpath("//tr[td]"):
        c=[" ".join(td.text_content().split()) for td in tr.xpath("td")]
        if len(c)>9:
            try:
                when=datetime.strptime(c[1],"%y%m%d %H:%M:%S.%f")
            except ValueError:
                continue
            events.append((when,c[2],c[8],c[0],c[3]))
    return sorted(events)


def compare(rows,gbm):
    epoch=Time(56658,0.000777592592592593,format="mjd",scale="tt")
    # TimeDelta uses SI seconds through leap offsets. MET already includes
    # event TIMEZERO. Subtracting the two UTC instants compares local clocks.
    fermi=(Time([entry[0] for entry in gbm],scale="utc").tt-epoch).sec
    matches=[]
    for r in rows:
        base=float(r["predicted_met_tt_s"])
        for shift in SHIFTS:
            t=base+shift
            lo,hi=np.searchsorted(fermi,[t-1,t+1],side="left")
            for i in range(lo,hi):
                utc,trig,source,met,duration=gbm[i]
                matches.append(dict(obsid=r["obsid"],source=r["source"],number=r["number"],
                   shift_s=shift,gbm_met=met,gbm_event_utc=utc.isoformat(timespec="milliseconds"),
                   gbm_trigger=trig,gbm_source=source,gbm_dur_ms=duration,
                   residual_s=round(fermi[i]-t,6),same_source=int(source==ALIASES.get(r["source"],r["source"]))))
    return matches


def summary(matches,radius):
    print(f"Within {radius} s (unique GBM detections, unique NICER rows, pairs, ObsIDs):")
    for shift in SHIFTS:
        m=[r for r in matches if r["shift_s"]==shift and abs(r["residual_s"])<=radius]
        s=[r for r in m if r["same_source"]]
        def count(group):
            return (len({r["gbm_met"] for r in group}),len({(r["obsid"],r["number"]) for r in group}),
                    len(group),len({r["obsid"] for r in group}))
        by_source=Counter(s0 for s0,_ in {(r["source"],r["gbm_met"]) for r in s})
        print(shift,"all",count(m),"same source",count(s),dict(sorted(by_source.items())))
    plus=[r for r in matches if r["shift_s"]==69.184 and r["same_source"] and abs(r["residual_s"])<=radius]
    print("Positive shift, other than ObsID 3020560101:",
          len({r["gbm_met"] for r in plus if r["obsid"]!="3020560101"}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("windows",type=Path)
    parser.add_argument("gbm_catalog",type=Path)
    parser.add_argument("--write",type=Path)
    args=parser.parse_args()
    with args.windows.open() as f:
        rows=list(csv.DictReader(f))
    events=read_gbm(args.gbm_catalog)
    if len(rows)!=1132 or len(events)!=1254:
        raise ValueError("Unexpected catalogue size")
    matches=compare(rows,events)
    for radius in (0.1,0.25,1):
        summary(matches,radius)
    if args.write:
        with args.write.open("w") as f:
            w=csv.DictWriter(f,list(matches[0]),lineterminator="\n")
            w.writeheader()
            w.writerows(matches)
