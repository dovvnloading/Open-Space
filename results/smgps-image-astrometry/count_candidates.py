"""Verify duplicate rows using both the source manuscript and CDS Table A.4."""

import csv
import hashlib
import io
import json
import tarfile
from urllib.request import urlopen

from harvest import HERE, get

TABLE = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2BA/693/A247/meerkatc&-out=recno,Name,GLON,GLAT,Rad,Rel&-out.max=unlimited"
MANUSCRIPT = "https://arxiv.org/src/2409.16607v1"


def main():
    groups = {}
    for line in get(TABLE).decode().splitlines():
        fields = [x.strip() for x in line.split("\t")]
        if len(fields) == 6 and fields[0].isdigit():
            groups.setdefault(tuple(fields[1:]), []).append(int(fields[0]))
    with urlopen(MANUSCRIPT, timeout=60) as response:
        data = response.read()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        table = tar.extractfile("MeerKAT_candidates.tex").read()
    repeats = [
        {"record_numbers": ids, "fields": list(fields)}
        for fields, ids in groups.items()
        if len(ids) > 1
    ]
    with (HERE / "headers.csv").open() as src:
        names = {
            row["name"]
            for row in csv.DictReader(src)
            if row["filename"].startswith("SMGPS_Candidate_")
        }
    summary = {
        "vizier_rows": sum(map(len, groups.values())),
        "distinct_full_records": len(groups),
        "distinct_table_names": len({x[0] for x in groups}),
        "distinct_image_names": len(names),
        "table_names_missing_images": sorted({x[0] for x in groups} - names),
        "image_names_missing_rows": sorted(names - {x[0] for x in groups}),
        "duplicate": repeats,
        "manuscript_table_sha256": hashlib.sha256(table).hexdigest(),
        "manuscript_occurrences": sum(
            line.startswith(b"G354.353$-$00.018") for line in table.splitlines()
        ),
    }
    (HERE / "candidate-count.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(summary)


if __name__ == "__main__":
    main()
