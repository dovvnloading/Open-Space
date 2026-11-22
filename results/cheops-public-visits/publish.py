"""Publish a bounded archive access result from the frozen evidence."""
import hashlib
import json
import re
from pathlib import Path
from collections import Counter

STUDY = Path(__file__).resolve().parent
OUT = STUDY/"outputs"


def catalog(sentence, label):
    path = Path("README.md")
    text = path.read_text()
    relative = f"results/{STUDY.name}/README.md"
    lines = text.splitlines()
    lines = [s for s in lines if not (s.startswith("- ") and relative in s)]
    section = lines.index("### Planets around other stars")
    lines[section+1:section+1] = ["",f"- {sentence} *{label}* · [Details]({relative})"]
    counts = Counter()
    for line in lines:
        if line.startswith("- "):
            tags = re.findall(r"(?<!\*)\*(Measured|Candidate|Prediction|Correction|Null result)\*(?!\*)",line)
            if tags:
                counts[tags[-1]]+=1
    line = lines[8]
    for pattern,count in [(r"(?<=has )\d+(?= entries)",sum(counts.values())),
                          (r"\d+(?= new measurements)",counts["Measured"]),
                          (r"\d+(?= candidate discoveries)",counts["Candidate"]),
                          (r"\d+(?= forecasts)",counts["Prediction"]),
                          (r"\d+(?= corrections)",counts["Correction"]),
                          (r"\d+(?= carefully recorded null)",counts["Null result"])]:
        line,n = re.subn(pattern,str(count),line,count=1)
        if n!=1:
            raise ValueError(f"Summary sentence format changed at {pattern}")
    lines[8] = line
    path.write_text("\n".join(lines)+"\n")


def checksums():
    lines = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name!="SHA256SUMS":
            lines.append(f"{hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}  {p.relative_to(OUT)}")
    (OUT/"SHA256SUMS").write_text("\n".join(lines)+"\n")


def publish():
    status = json.loads((OUT/"archive_status.json").read_text())
    probes = [json.loads(l) for l in (OUT/"archive_probes.jsonl").read_text().splitlines()]
    ready = status["verified_lightcurve"]
    count = status["visits"]
    tested = len(status["tested_file_keys"])
    if ready:
        claim = "The public archive returned a CHEOPS light curve"
        summary = (f"The inventory contains {count:,} visit identifiers. At least one anonymous download "
                   "yielded a FITS table with time and flux samples. The access run did not scan for transits.")
    else:
        claim = "The anonymous CHEOPS access test returned no verified light curves"
        summary = (f"The archive queries yielded {count:,} visit identifiers. Product requests covered "
                   f"{tested} distinct visits. Zero light curves reached the transit search. "
                   "No stars, epochs, signal depths or durations can be excluded. Injection recovery was unmeasured.")
        catalog("**The anonymous CHEOPS archive check returned no usable light curve.** "
                "The transit search stopped at data access. It sets no limit on missed planets or eclipses.",
                "Null result")
    entries=[]
    for p in probes:
        if "url" not in p and not p.get("error"):
            continue
        state = str(p.get("http",""))+" "+str(p.get("error",""))
        state = state.replace("|","/").replace("\n"," ")[:220]
        value=(p["operation"],state.strip())
        if value not in entries:
            entries.append(value)
    rows = "\n".join("| "+a.replace("|","/")+" | "+b+" |" for a,b in entries)
    report = f"""# {claim}

{summary}

## What the archive returned

| Request | HTTP status or recorded failure |
| --- | --- |
{rows}

The inventory records which service supplied each visit. The probes record exact queries, response status and SHA-256 digests of retrieved inputs. A successful metadata or portal request alone cannot establish photometry access. The recorded status was `{status['status']}`. Credential challenges, HTTP service failures and payload parsing failures appear separately in the receipts. The conclusion applies to these requests at the time of the run.

[Archive status](outputs/archive_status.json) · [Visit inventory](outputs/inventory.jsonl) · [Request and input receipts](outputs/archive_probes.jsonl) · [Checksums](outputs/SHA256SUMS)

## Reproduction

In the study folder, the commands used by the archive job were

```sh
python -m pip install -r requirements.txt
python archive_io.py inventory outputs
```

The archive is live and can change. Verify the deposited snapshot offline from the study folder with

```sh
(cd outputs && sha256sum -c SHA256SUMS)
```
"""
    (STUDY/"README.md").write_text(report)
    checksums()


if __name__=="__main__":
    publish()
