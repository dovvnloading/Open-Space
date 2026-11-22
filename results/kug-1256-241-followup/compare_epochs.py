"""Summarize nearest same-band observing seasons around the slew and eRASS1."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

OUT = Path(__file__).resolve().parent / "outputs"

def main():
    with (OUT/"timeline.csv").open() as f:
        rows = list(csv.DictReader(f))
    data=defaultdict(list)
    epochs={}
    for r in rows:
        if r["archive"] in ("XMM slew","eROSITA eRASS1") and r["mjd"]:
            epochs[r["archive"]]=(float(r["mjd"])+float(r["end_mjd"] or r["mjd"]))/2
        if r["kind"]=="detection" and r["mjd"] and r["value"] and r["unit"].startswith("mag"):
            data[(r["archive"],r["band"],r["unit"])].append((float(r["mjd"]),float(r["value"])))
    output=[]
    for (a,b,u),series in sorted(data.items()):
        for label,day in sorted(epochs.items()):
            # Symmetric windows maintain the same cross-band timing tolerance.
            near=[(t,m) for t,m in series if abs(t-day)<=45]
            if near:
                v=median(m for t,m in near)
                output.append([label,round(day,5),a,b,u,len(near),
                               round(min(t for t,m in near),5),round(max(t for t,m in near),5),
                               round(v,4),round(median(abs(m-v) for t,m in near),4)])
    with (OUT/"optical-ir-near-xray.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["xray","xray_mjd","archive","band","unit","n_within_45d","first_mjd","last_mjd","median","mad"])
        w.writerows(output)


if __name__=="__main__":
    main()
