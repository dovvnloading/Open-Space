"""Check the archive replay against the frozen epoch-by-epoch measurements."""
import csv
import json
from pathlib import Path
import sys


def main():
    root = Path(sys.argv[1])
    data = json.loads((root / "measurements.json").read_text())
    with Path(__file__).with_name("observations.csv").open() as f:
        rows = list(csv.DictReader(f))
    names = {"EROS red":"EROS R","EROS blue":"EROS B","EROS red before OGLE III":"EROS before OGLE-III",
             "OGLE III I":"OGLE-III I","OGLE IV I":"OGLE-IV I","OGLE IV V":"OGLE-IV V"}
    for row in rows:
        k,label = row["target"],row["measurement"]
        if label.startswith("OGLE I block"):
            result = data["shortening"][k]["days"][int(label.split()[-1])-1]
        else:
            result = data["dominance" if k=="046289" else "shortening"][k][names.get(label,label)]
        assert result["days"] == int(row["occupied_days"]), row
        assert abs(float(row["period_days"])-result["period_d"]) < .02, row
        assert abs(float(row["variance_fraction"])-result["variance_fraction"]) < .003, row
    screen = root / "screen.csv"
    if screen.exists():
        with screen.open() as f:
            ledger = list(csv.DictReader(f))
        assert len(ledger)==52 and sum(r["fit_III"]!="404" and r["fit_IV"]!="404" for r in ledger)==45
    print(f"Verified {len(rows)} independently fitted period records.")


if __name__ == "__main__":
    main()
