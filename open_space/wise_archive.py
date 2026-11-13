"""Bounded, replayable cone queries against two public IRSA single-source catalogs."""

import csv
import hashlib
import io
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .archive import write_immutable

IRSA = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
COLUMNS = {
    "neowise": (
        "source_id",
        "scan_id",
        "frame_num",
        "mjd",
        "ra",
        "dec",
        "sigra",
        "sigdec",
        "w1mpro",
        "w1sigmpro",
        "w1snr",
        "w1rchi2",
        "w1sat",
        "cc_flags",
        "qual_frame",
        "qi_fact",
        "moon_masked",
        "saa_sep",
        "sso_flg",
        "nb",
        "na",
        "allwise_cntr",
    ),
    "gaia": ("source_id", "ra", "dec", "pmra", "pmdec", "ref_epoch", "ruwe"),
}
TABLES = {"neowise": "neowiser_p1bs_psd", "gaia": "gaia_dr3_source"}
# IRSA applies these before transmitting the rows; all other cuts are replayed locally.
WISE_CONSTRAINT = "w1snr >= 8 and qual_frame > 0"
MAX_ROWS = 60000
MAX_BYTES = 50_000_000


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def field_params(field):
    ra, dec, radius = (float(field[k]) for k in ("ra_deg", "dec_deg", "radius_arcmin"))
    if not (0 <= ra < 360 and -85 <= dec <= 85 and 0 < radius <= 10):
        raise ValueError("Use finite ICRS coordinates, |Dec| <= 85 deg, and 0 < radius <= 10 arcmin")
    if field.get("role") not in ("exploration", "positive_control") or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_-]{0,64}", field["id"]
    ):
        raise ValueError("Invalid WISE field ID or role")
    return ra, dec, radius


def query_url(field, kind):
    ra, dec, radius = field_params(field)
    # The 15 arcsec Gaia margin admits tracks drifting across the WISE boundary.
    args = {
        "catalog": TABLES[kind],
        "spatial": "Cone",
        "objstr": f"{ra:.8f} {dec:+.8f} eq",
        "radius": f"{radius * 60 + (15 if kind == 'gaia' else 0):.4f}",
        "radunits": "arcsec",
        "outfmt": 1,
        "selcols": ",".join(COLUMNS[kind]),
        "outrows": MAX_ROWS + 1,
    }
    if kind == "neowise":
        args["constraints"] = WISE_CONSTRAINT
    return IRSA + "?" + urlencode(args)


def canonical_table(payload, kind):
    """Keep selected text fields without numeric precision loss or transient HTTP metadata."""
    lines = payload.decode("ascii").splitlines()
    try:
        num = int(next(s.split("=", 1)[1] for s in lines if s.startswith("\\RowsRetrieved")))
        header_idx = next(i for i, line in enumerate(lines) if line.startswith("|"))
    except (ValueError, StopIteration) as exc:
        raise ValueError(f"IRSA returned an unexpected table: {payload[:220]!r}") from exc
    if num > MAX_ROWS:
        raise ValueError("IRSA field exceeds the row limit; split it into smaller fields")
    heading = lines[header_idx]
    marks = [i for i, c in enumerate(heading) if c == "|"]
    bounds = {heading[a + 1 : b].strip(): (a + 1, b) for a, b in zip(marks, marks[1:])}
    if not set(COLUMNS[kind]).issubset(bounds):
        raise ValueError("IRSA result is missing required columns")
    rows = []
    for line in lines[header_idx + 1 :]:
        if not line.strip() or line.startswith(("|", "\\")):
            continue
        rows.append(tuple(line[bounds[c][0] : bounds[c][1]].strip() for c in COLUMNS[kind]))
    if len(rows) != num:
        raise ValueError(f"Incomplete IRSA result ({len(rows)} of {num} rows)")
    rows.sort()
    if len({v[0] for v in rows}) != len(rows):
        raise ValueError("IRSA returned duplicate detection or source identifiers")
    out = io.StringIO(newline="")
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(COLUMNS[kind])
    writer.writerows(rows)
    return out.getvalue().encode("ascii")


def acquire(field, cache_dir, kind):
    url = query_url(field, kind)
    root = Path(cache_dir) / field["id"] / digest(url.encode())[:16]
    root.mkdir(parents=True, exist_ok=True)
    path = root / (kind + ".csv")
    raw_path = root / (kind + ".tbl")
    expected = field.get("input_sha256", {}).get(kind)
    receipt = root / (kind + ".json")
    if not receipt.exists():
        if not raw_path.exists():
            for attempt in range(3):
                try:
                    with urlopen(
                        Request(url, headers={"User-Agent": "Open-Space-Archive-Research/1.0"}),
                        timeout=150,
                    ) as response:
                        data = response.read(MAX_BYTES + 1)
                    if len(data) > MAX_BYTES:
                        raise ValueError("IRSA response exceeded 50 MB")
                    normalized = canonical_table(data, kind)
                    if expected and digest(normalized) != expected:
                        raise ValueError("IRSA normalized input has changed since the published run")
                    write_immutable(raw_path, data)
                    break
                except (HTTPError, URLError, TimeoutError) as exc:
                    if (
                        isinstance(exc, HTTPError) and exc.code < 500 and exc.code != 429
                    ) or attempt == 2:
                        raise
                    time.sleep(1 + 2 * attempt)
        data = raw_path.read_bytes()
        normalized = canonical_table(data, kind)
        if expected and digest(normalized) != expected:
            raise ValueError("Cached IRSA input differs from the requested version")
        write_immutable(path, normalized)
        meta = {
            "archive": "NASA/IPAC IRSA",
            "table": TABLES[kind],
            "url": url,
            "sha256": digest(normalized),
            "response_sha256": digest(data),
            "row_count": len(normalized.splitlines()) - 1,
        }
        # HTTP headers embedded by Gator include retrieval times and transient paths.
        write_immutable(receipt, (json.dumps(meta, sort_keys=True, indent=2) + "\n").encode())
    meta = json.loads(receipt.read_bytes())
    data = path.read_bytes()
    if (
        meta["url"] != url
        or digest(data) != meta["sha256"]
        or (expected and digest(data) != expected)
        or digest(raw_path.read_bytes()) != meta["response_sha256"]
    ):
        raise ValueError("IRSA cached input or acquisition metadata failed verification")
    # Run evidence references the canonical science table; the raw receipt stays in the cache.
    stable = {k: meta[k] for k in ("archive", "table", "url", "sha256", "row_count")}
    return list(csv.DictReader(io.StringIO(data.decode("ascii")))), stable
