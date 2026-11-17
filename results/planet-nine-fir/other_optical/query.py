"""NOIRLab current DECam archive searches near the moving Chen hypotheses."""

from pathlib import Path
import hashlib, json, sys
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import datetime
from concurrent.futures import ThreadPoolExecutor

HERE = Path(__file__).resolve().parent
CACHE = HERE / "_cache"
CACHE.mkdir(exist_ok=True)


def acquire(name, url, spec=None):
    rec = dict(url=url, utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    req = Request(url)
    if spec is not None:
        rec["payload"] = spec
        req = Request(
            url,
            data=json.dumps(spec).encode(),
            headers={"Content-Type": "application/json"},
        )
    try:
        with urlopen(req, timeout=100) as res:
            data = res.read()
            (CACHE / (name + ".raw")).write_bytes(data)
            rec.update(
                status=res.status,
                sha256=hashlib.sha256(data).hexdigest(),
                size=len(data),
            )
    except Exception as e:
        rec["error"] = repr(e)
    (CACHE / (name + ".json")).write_text(json.dumps(rec, indent=2) + "\n")
    print(name, rec.get("size"), rec.get("error"), flush=True)


def field_query(item):
    name, a, d = item
    q = {
        "outfields": [
            "archive_filename",
            "md5sum",
            "dateobs_center",
            "caldat",
            "ra_center",
            "dec_center",
            "ifilter",
            "exposure",
            "proposal",
            "release_date",
            "OBJECT",
        ],
        "search": [
            ["ra_center", a - 3, a + 3],
            ["dec_center", d - 2.9, d + 2.9],
            ["caldat", "2012-01-01", "2026-11-17"],
            ["instrument", "decam"],
            ["proc_type", "instcal"],
            ["prod_type", "image"],
        ],
    }
    acquire(
        name,
        "https://astroarchive.noirlab.edu/api/adv_search/find/?rectype=file&limit=20000",
        q,
    )


if __name__ == "__main__":
    with ThreadPoolExecutor(2) as pool:
        list(
            pool.map(
                field_query, [("Chen1", 42.676, -15.021), ("Chen2", 45.297, -16.711)]
            )
        )
