"""Convert the thirty unchanged IRSA source responses to one smaller reproducible input."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


def build(pairs, directory, destination):
    fields = ["source_id", "oid", "hjd", "mag", "magerr", "catflags", "filtercode", "ra", "dec"]
    manifest = []
    with destination.open("wb") as file, gzip.GzipFile(filename="",fileobj=file,mode="wb",mtime=0) as archive:
        with io.TextIOWrapper(archive,encoding="utf-8",newline="") as sink:
            writer = csv.DictWriter(sink, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for pair in csv.DictReader(pairs.open()):
                for i in (1,2):
                    sid = pair[f"source_id{i}"]
                    raw = (directory / (sid+".csv")).read_bytes()
                    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
                    for row in rows:
                        writer.writerow({"source_id":sid, **{field:row[field] for field in fields[1:]}})
                    manifest.append(dict(source_id=sid, rows=len(rows), sha256=hashlib.sha256(raw).hexdigest()))
    return manifest


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("pairs","directory","destination","manifest"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    a.manifest.write_text(json.dumps(build(a.pairs,a.directory,a.destination),indent=2)+"\n")
