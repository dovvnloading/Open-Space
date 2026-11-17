"""Test the blind W2 links on NEOWISE exposures after the last unTimely epoch."""
import csv
import json
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"widebatch"))
import screen

def main():
    screen.BASE=ROOT/"l1b"
    for end in ["cache","retained"]:(screen.BASE/end).mkdir(exist_ok=True,parents=True)
    rows=list(csv.DictReader((ROOT/"tracks.csv").open()))
    def run(row):
        r=screen.one(row,heldout_mjd=59215.)
        r["role"]="known Y in control field" if row["source_name"]=="track1342m076_1" else "control-field link" if row["tile"]=="1342m076" else "blind link"
        return r
    with ThreadPoolExecutor(4) as ex:results=list(ex.map(run,rows))
    with (screen.BASE/"summary.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
    manifest=[]
    for row in rows:
        m=json.loads((screen.BASE/"retained"/(row["source_name"]+".json")).read_text())
        manifest.append({k:v for k,v in m.items() if k in ("source","archive_url","raw_ipac_sha256","raw_science_rows_sha256","retained_csv_sha256","field_mjd","heldout_mjd","heldout_2arcsec","quality_selected_epochs","prior_selected_five_parameter_fit")})
    (screen.BASE/"manifest.json").write_text(json.dumps(manifest,indent=1)+"\n")
if __name__=="__main__":main()
