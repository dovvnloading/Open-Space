"""Guest TNS cone queries without an API key; make access failures explicit.

Public search page parameters follow https://www.wis-tns.org/content/faq
and its downloadable per-page CSV query convention. This test does NOT
establish novelty when the server blocks, redirects, or returns HTML.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import csv
import io
import json

POSITIONS = {"COSMOS-022": (150.122033, 2.323386208),
             "GOODS-S 185": (53.188672, -27.812554)}


def check(label, ra, dec, view):
    args = {"ra": str(ra), "decl": "%+.9f" % dec,
            "radius": "10", "coords_unit": "arcsec",
            "public": "all", "isTNS_AT": "all", "format": view}
    if view == "csv":
        args.update({"num_page": "100", "page": "0"})
    url = "https://www.wis-tns.org/search?" + urlencode(args)
    record = {"candidate": label, "view": view, "url": url,
              "cone_radius_arcsec": 10}
    try:
        with urlopen(Request(url, headers={
            "User-Agent": "Mozilla/5.0 (compatible; independent public-catalog review)"}), timeout=25) as r:
            body = r.read(2_000_000)
            record.update(status=r.status, final_url=r.url,
                          content_type=r.headers.get("Content-Type"),
                          response_bytes_read=len(body),
                          response_sha256=sha256(body).hexdigest())
            if view == "csv" and "csv" in (r.headers.get("Content-Type") or "").lower():
                rows = list(csv.reader(io.StringIO(body.decode("utf-8", errors="replace"))))
                record["returned_csv_rows"] = max(0, len(rows) - 1)
                record["outcome"] = "CSV retrieved; inspect positions independently, not proof of novelty"
            else:
                record["outcome"] = "Search page retrieved, not machine-verified for counterpart count"
    except HTTPError as e:
        record.update(status=e.code, outcome="TNS HTTP access denied/unavailable",
                      response_excerpt=e.read(200).decode("utf-8", errors="replace"))
    except (URLError, TimeoutError, OSError) as e:
        record.update(status=None, outcome="Network/TNS request failed",
                      error=str(e)[:220])
    return record


if __name__ == "__main__":
    jobs = [(label, ra, dec, view) for label, (ra, dec) in POSITIONS.items()
            for view in ("html", "csv")]
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(lambda args: check(*args), jobs))
    print(json.dumps({"retrieved_utc": datetime.now(timezone.utc).isoformat(),
                      "source": "TNS public guest search HTML and FAQ page-CSV route",
                      "results": records,
                      "caution": "HTTP errors and inaccessible CSV do not establish absence of a TNS entry; previously frozen brokers still have incomplete coverage"},
                     indent=2))
