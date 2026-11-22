"""Gather literature and counterpart evidence for the survey's passing sources."""
import json
import sys
from pathlib import Path

import vet

work = Path(sys.argv[1])
leads = json.loads((work / "leads.json").read_text())
for lead in leads:
    lead["aliases"] = [str(lead[k]) for k in ("assoc",) if lead.get(k)]
destination = work / "review"
destination.mkdir(parents=True, exist_ok=True)
result = vet.collect(leads, destination)
(work / "review-summary.json").write_text(json.dumps(result, indent=2) + "\n")
