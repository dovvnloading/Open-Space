"""Fetch the pinned archival images and query outputs listed in inputs.json."""

import hashlib, json, time
from pathlib import Path
from urllib.request import Request, urlopen
from concurrent.futures import ThreadPoolExecutor

here = Path(__file__).resolve().parent
cache = here / "_cache"
cache.mkdir(exist_ok=True)


def fetch(row):
    path = cache / (row["label"] + ".raw")
    receipt = cache / (row["label"] + ".json")
    if path.exists():
        with path.open("rb") as f:
            if hashlib.file_digest(f, "sha256").hexdigest() == row["sha256"]:
                receipt.write_text(json.dumps(row))
                return
    for retry in range(3):
        try:
            req = Request(row["url"])
            if "payload" in row:
                req = Request(
                    row["url"],
                    data=json.dumps(row["payload"]).encode(),
                    headers={"Content-Type": "application/json"},
                )
            with urlopen(req, timeout=150) as res, path.with_suffix(".part").open(
                "wb"
            ) as out:
                digest = hashlib.sha256()
                while data := res.read(1 << 20):
                    out.write(data)
                    digest.update(data)
            if digest.hexdigest() != row["sha256"]:
                raise ValueError(
                    "Different bytes for "
                    + row["label"]
                    + "; inspect the query and data revision before using."
                )
            path.with_suffix(".part").replace(path)
            receipt.write_text(json.dumps(row))
            return
        except (OSError, TimeoutError):
            if retry == 2:
                raise
            time.sleep(2 * (retry + 1))


if __name__ == "__main__":
    with ThreadPoolExecutor(3) as e:
        list(
            e.map(
                fetch,
                (
                    r
                    for r in json.loads((here / "inputs.json").read_text())
                    if r["required_for_reduction"]
                ),
            )
        )
