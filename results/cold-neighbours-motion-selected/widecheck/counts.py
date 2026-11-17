"""Summarize the six-visit screen without treating positional matches as IDs."""
import csv,json,statistics
from collections import Counter
from pathlib import Path
BASE=Path(__file__).resolve().parent
def main():
    with (BASE/"summary.csv").open() as f:r=list(csv.DictReader(f))
    open_alerts=[s for s in r if s["priority"]=="1" and s["listed_match"]=="0"]
    lit=[s for s in r if s["listed_match"]=="1"]
    m=json.loads((BASE/"manifest.json").read_text())
    epochs=sum(int(x["length"]) for x in m["catalogues"])
    d=dict(sources=len(r),tiles=len({x["source_id"][:8] for x in r}),files=len(m["catalogues"]),compressed_catalogue_bytes=epochs,
      files_per_source=dict(Counter(x["n_files"] for x in r)),usable_neighbor_epochs=dict(Counter(x["local_depth_epochs"] for x in r)),
      local_depth_at_least_three=sum(int(x["depth_reaches_w2"])>=3 for x in r),
      fixed_window_counts=dict(Counter(x["moving"] for x in r)),
      passing=sum(s["priority"]=="1" for s in r),
      literature_alerts=len(lit),passing_with_literature_alert=sum(s["priority"]=="1" for s in lit),
      passing_without_literature_alert=len(open_alerts),unmatched_passes=[s["source_name"] for s in open_alerts],
      rotated_passes={k:sum(x[k+"_passes"]=="1" for x in r) for k in ["rot90","rot180","rot270"]})
    (BASE/"counts.json").write_text(json.dumps(d,indent=2,sort_keys=True)+"\n");print(d)
if __name__=="__main__":main()
