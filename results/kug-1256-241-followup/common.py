"""One target, UTC/MJD timeline and per-request provenance."""
import csv
import hashlib
import json
import math
import threading
from datetime import datetime, timezone
from pathlib import Path

import requests
from astropy.time import Time

RA, DEC = 194.714898, 23.923199
HOST_RA, HOST_DEC = 194.714542, 23.924111
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
PREVIOUS = Path("results/xmmsl2-erass1-rass")
OUT.mkdir(parents=True, exist_ok=True)
lock = threading.Lock()
events, coverage, receipts = [], [], []


def fetch(key, url, params=None, *, data=None, timeout=45):
    """Hash the actual uncompressed response bytes that each parser consumes."""
    item = dict(key=key, url=url, params=params, method="POST" if data is not None else "GET",
                data=data, utc=datetime.now(timezone.utc).isoformat())
    try:
        response = requests.request(item["method"], url, params=params, data=data,
                                    headers={"User-Agent": "Open-Space/1.0 (public astronomy research; github.com/dovvnloading/Open-Space)"},
                                    timeout=timeout)
        item.update(http=response.status_code, final_url=response.url, bytes=len(response.content),
                    sha256=hashlib.sha256(response.content).hexdigest())
        response.raise_for_status()
        return response.content
    except Exception as exc:
        item["error"] = str(exc)[:600]
        raise
    finally:
        with lock:
            receipts.append(item)


def input_file(path):
    path = Path(path)
    size = path.stat().st_size
    with path.open("rb") as file:
        digest = hashlib.file_digest(file, "sha256").hexdigest()
    with lock:
        receipts.append(dict(key=path.name, method="git input", path=str(path),
                             sha256=digest, bytes=size))
    return path


def event(archive, band, mjd, value="", error="", unit="", kind="detection",
          source="", note="", end_mjd=""):
    def time_pair(value):
        if value in ("", None):
            return "", ""
        day = float(value)
        if not math.isfinite(day):
            return "", ""
        return day, Time(day, format="mjd", scale="utc").isot
    start, utc_start = time_pair(mjd)
    end, utc_end = time_pair(end_mjd)
    entry = dict(mjd=start, utc=utc_start, end_mjd=end, end_utc=utc_end,
                 archive=archive, band=band, kind=kind, value=value, error=error,
                 unit=unit, source=source, note=note)
    with lock:
        events.append(entry)


def status(archive, state, detail):
    with lock:
        coverage.append(dict(archive=archive, state=state, detail=str(detail)))


def write_csv(path, rows, columns):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def save():
    timeline = sorted(events, key=lambda r: (r["mjd"] == "", r["mjd"] or 0,
                                              str(r["archive"]), str(r["band"]), str(r["source"])))
    write_csv(OUT / "timeline.csv", timeline,
              "mjd utc end_mjd end_utc archive band kind value error unit source note".split())
    write_csv(OUT / "coverage.csv", sorted(coverage, key=lambda r: (r["archive"], r["state"])),
              ["archive", "state", "detail"])
    (OUT / "inputs.jsonl").write_text("".join(json.dumps(r, default=str, ensure_ascii=False) + "\n"
                                             for r in sorted(receipts, key=lambda r: str(r["key"]))))
