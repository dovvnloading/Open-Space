"""Run the held-out W2 count criterion at all wide literature-alert positions; compare by W2 brightness."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import sys
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"widebatch"));import screen

def main():
    p=argparse.ArgumentParser();p.add_argument("--selected",type=Path,required=True);p.add_argument("--matches",type=Path,required=True);a=p.parse_args()
    (ROOT/"cache").mkdir(exist_ok=True);(ROOT/"retained").mkdir(exist_ok=True)
    s=pd.read_csv(a.selected);matches=pd.read_csv(a.matches)
    names=set(matches.source_name);w=s[s.source_name.isin(names)].copy()
    screen.BASE=ROOT
    def run(row):
        rec=screen.one(row)
        subset=matches[matches.source_name==row["source_name"]]
        rec["catalogue_alerts"]="|".join(sorted(set(subset.catalogue)))
        rec["aliases"]="|".join(sorted(set(subset.matched_label)))
        rec["role"]="literature positional alert"
        return rec
    with ThreadPoolExecutor(2) as ex:results=list(ex.map(run,w.to_dict('records')))
    summary=pd.DataFrame(results)
    summary.to_csv(ROOT/"summary.csv",index=False)
    marks=summary.triage.str.startswith("track-supported")
    grouped=[]
    for low,high in [(0.,14.),(14.,14.5),(14.5,15.),(15.,15.501)]:
        cell=(summary.w2>=low)&(summary.w2<high)
        grouped.append(dict(w2_min=low,w2_max_exclusive=high,n_alerts=int(cell.sum()),n_supported=int((marks&cell).sum())))
    pd.DataFrame(grouped).to_csv(ROOT/"recovery.csv",index=False)
    # Frozen archive queries and source-row digests; all centroids can be recreated from archive endpoints.
    fields=["source","archive_url","raw_ipac_sha256","raw_science_rows_sha256","retained_csv_sha256","field_mjd","heldout_mjd","heldout_2arcsec"]
    records=[]
    for row in results:
        data=json.loads((ROOT/"retained"/(row["source"]+".json")).read_text())
        records.append({key:data[key] for key in fields})
    (ROOT/"inputs.json").write_text(json.dumps(records,indent=1)+"\n")
    print(grouped)
if __name__=="__main__":main()
