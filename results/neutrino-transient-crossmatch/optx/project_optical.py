"""Project the complete frozen inputs from receipts.json; arg is cache directory with receipt keys as filenames. Requires fastparquet and pandas."""
import csv
import sys
from pathlib import Path
import fastparquet
import pandas as pd
from optical import DIR

cache = Path(sys.argv[1])
columns = ["objid","name_prefix","name","ra","declination","type","discoverydate","reporting_group","source_group","internal_names","lastmodified"]
frames = [fastparquet.ParquetFile(cache / f"p{i}").to_pandas(columns=columns) for i in range(12)]
d = pd.concat(frames, ignore_index=True)
assert len(d) == 213398 and d.objid.nunique() == len(d)
sel = d.type.str.fullmatch("TDE(-H|-He|-H-He|-featureless)?", na=False)
d.loc[sel & (d.discoverydate >= "2024") & (d.discoverydate <= "2026-11-17 12:09:02"), columns].sort_values("objid").to_csv(DIR/"data"/"tns_selected.csv",index=False)
for key, dest in [("tdes","tdecat"),("btscsv","bts"),("lionscsv","lions")]:
    data = list(csv.DictReader((cache/key).read_text(encoding="utf-8-sig").strip().splitlines()))
    if key=="tdes":
        rows=[x for x in data if x["Discovery date (UT)"] >= "2024"]
    elif key=="btscsv":
        rows=[x for x in data if x["discdate"] >= "2024" and x["type"].startswith("TDE")]
    else:
        rows=[x for x in data if x["classification"]=="Tidal Disruption Event"]
    with (DIR/"data"/f"{dest}_selected.csv").open("w",newline="") as out:
        w=csv.DictWriter(out,data[0],lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
