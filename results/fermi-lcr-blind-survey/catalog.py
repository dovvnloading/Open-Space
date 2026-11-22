"""Add one survey entry to the current catalog and count the actual labels."""
import json
import re
import sys
from collections import Counter
from pathlib import Path

path = Path("README.md")
lines = path.read_text().splitlines(keepends=True)
entry = json.loads(Path(sys.argv[1]).read_text())
topic = "### Galaxies and the distant universe"
link = "results/fermi-lcr-blind-survey/README.md"
bullet = f"- **{entry['title']}** {entry['detail']} *{entry['label']}* · [Details]({link})\n"
previous = [i for i, line in enumerate(lines) if link in line and line.startswith("- ")]
if len(previous) == 1:
    lines[previous[0]] = bullet
elif not previous:
    at = next(i for i,line in enumerate(lines) if line.rstrip() == topic)
    start = next(i for i in range(at+1, len(lines)) if lines[i].startswith("- "))
    lines.insert(start, bullet)
else:
    raise ValueError("Duplicate Fermi catalog entries")
labels = Counter()
for line in lines:
    if re.match(r"^\s*[-*+]\s+", line):
        matches = re.findall(r"(?<!\*)\*(Measured|Candidate|Prediction|Correction|Null result)\*(?!\*)", line)
        if matches:
            labels[matches[-1]] += 1
line = lines[8]
for pattern, n in [(r"(?<=has )\d[\d,]*(?= entries)", sum(labels.values())),
                   (r"\d[\d,]*(?= new measurements)", labels["Measured"]),
                   (r"\d[\d,]*(?= candidate discoveries)", labels["Candidate"]),
                   (r"\d[\d,]*(?= forecasts)", labels["Prediction"]),
                   (r"\d[\d,]*(?= corrections)", labels["Correction"]),
                   (r"\d[\d,]*(?= carefully recorded)", labels["Null result"])]:
    line, found = re.subn(pattern, str(n), line, count=1)
    if found != 1:
        raise ValueError("Unrecognized count wording on line 9")
lines[8] = line
path.write_text("".join(lines))
