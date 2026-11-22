"""Exercise live extraction on published BAT/MAXI sources before processing all rows."""
import json
from pathlib import Path
import sys

from ingest import load_curve

out = Path(sys.argv[1])
rows = json.loads((out / "inventory.json").read_text())["sources"]
(out / "pilot" / "curves").mkdir(exist_ok=True, parents=True)
control_names = {"crab", "herx-1", "smcx-1", "cenx-3", "velax-1"}
pilot = [r for r in rows if r["id"].lower() in control_names]
report = [{"id": r["id"], "name": r["name"], "maxi_id": r.get("maxi_id"),
           **{k:load_curve(r, k, out / "pilot") for k in ("bat", "maxi")}} for r in pilot]
(out / "controls-ingest.json").write_text(json.dumps(report, indent=2)+"\n")
for r in report:
    print(r["id"], [(k, r[k]["status"], r[k].get("rows_used"), r[k].get("first_mjd"), r[k].get("last_mjd"),r[k].get("attempts", [])[-2:]) for k in ("bat", "maxi")])
if len(pilot) < 2 or sum(all(r[k]["status"] == "ok" for k in ("bat", "maxi")) for r in report) < 2:
    raise SystemExit("Fewer than two control sources yielded both instruments.")
