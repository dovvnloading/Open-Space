"""Inventory the official monitor rosters and fetch their independently measured daily rates.

Only the workflow runner should invoke this module. Raw SHA-256 digests refer to
the HTTP entity after transfer decompression and before any FITS/ASCII decoding.
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import io
import json
import math
import re
import threading
import time
from pathlib import Path
from urllib.parse import unquote, urljoin, urlparse

import numpy as np
import requests
from astropy.io import fits
from bs4 import BeautifulSoup

BAT = "https://swift.gsfc.nasa.gov/results/transients/"
MAXI = "https://maxi.riken.jp/"
STATE = threading.local()


def fetch(url):
    if not hasattr(STATE, "session"):
        STATE.session = requests.Session()
        STATE.session.headers["User-Agent"] = "Open-Space BAT-MAXI survey (public light curves)"
    for attempt in range(3):
        try:
            r = STATE.session.get(url, timeout=(9, 42))
            r.raise_for_status()
            if not r.content:
                raise ValueError("empty response")
            return r.content, {"url": r.url, "sha256": hashlib.sha256(r.content).hexdigest(),
                               "bytes": len(r.content), "last_modified": r.headers.get("Last-Modified")}
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
            if (isinstance(exc, requests.HTTPError) and exc.response.status_code not in
                    (408, 429, 500, 502, 503, 504)) or attempt == 2:
                raise
            time.sleep(2 ** attempt)


def clean(value):
    return " ".join(str(value).split())


def key(name):
    # Preserve '+' and '-' since these distinguish many coordinate-based names.
    return re.sub(r"[^a-z0-9+\-]", "", name.casefold().replace("−", "-").replace("–", "-"))


def dedup(values):
    return list(dict.fromkeys(v for v in values if v))


def source_id(slug):
    return re.sub(r"[^A-Za-z0-9_.+\-]", "_", unquote(slug))


def page_links(data, base):
    soup = BeautifulSoup(data, "html.parser")
    for anchor in soup.find_all("a", href=True):
        url = urljoin(base, anchor["href"].strip()).split("#")[0]
        yield anchor, url, clean(anchor.get_text(" ", strip=True))


def row_cells(row):
    return [clean(x.get_text(" ", strip=True)) for x in row.find_all(["td", "th"], recursive=False)]


def degrees(ra, dec):
    try:
        xy = [float(ra.strip()), float(dec.strip().replace("−", "-"))]
    except ValueError:
        return None
    return xy if 0 <= xy[0] < 360 and -90 <= xy[1] <= 90 else None


def bat_index(data, url):
    rows = []
    soup = BeautifulSoup(data, "html.parser")
    # Every data row of the BAT root table is enumerated independently by the
    # archive. Reconcile against that enumeration to detect missing hrefs, a
    # changed table schema, accidental duplicate IDs, or a truncated response.
    indices = []
    trs = soup.find_all("tr")
    for tr in trs:
        cells = tr.find_all(["td", "th"], recursive=False)
        if not cells:
            continue
        number = clean(cells[0].get_text(" ", strip=True))
        if not number.isdigit():
            continue
        indices.append(int(number))
        if len(cells) < 6:
            raise ValueError(f"short BAT roster row {number}")
        anchor = cells[1].find("a", href=True)
        if anchor is None:
            raise ValueError(f"missing BAT product for row {number}")
        dest = urljoin(url, anchor["href"].strip()).removesuffix("/").removesuffix("/index.html")
        label = clean(anchor.get_text(" ", strip=True))
        slug = unquote(urlparse(dest).path.rsplit("/", 1)[-1])
        words = [clean(c.get_text(" ", strip=True)) for c in cells[2:6]]
        xy = degrees(*words[:2])
        if not xy:
            raise ValueError(f"invalid BAT position for row {number}")
        rows.append({"id": source_id(slug), "name": label, "bat_page": dest+"/",
                     "aliases": dedup([label, slug] + [a for a in re.split(r",\s*|\s+/\s+", words[2])
                                                        if a.strip() not in {"", "-", "--"}]),
                     "position": xy, "source_type": words[3], "archive_row": int(number)})
    if not indices or indices != list(range(1, len(indices)+1)) or len({b["id"] for b in rows}) != len(rows):
        raise ValueError("BAT enumeration is incomplete or contains nonunique product identifiers")
    return rows, []


def maxi_index(data, url):
    # The RA-sorted and name-sorted tables mirror each other. Use the latter and
    # collect all official products for each row, including older releases where
    # the standard product is absent. The four survey versions are not independent.
    found = {}
    soup = BeautifulSoup(data, "html.parser")
    for tr in soup.find_all("tr"):
        cells = tr.find_all("td", recursive=False)
        if len(cells) < 3:
            continue
        label, coord = [clean(c.get_text(" ", strip=True)) for c in cells[:2]]
        position = degrees(*coord.split(",")) if coord.count(",") == 1 else None
        if not position:
            continue
        for anchor in tr.find_all("a", href=True):
            dest = urljoin(url, anchor["href"].strip()).replace("http://maxi.riken.jp/", MAXI)
            m = re.search(r"/(J\d{4}[+\-]\d{3,4}[A-Za-z0-9]*)/", unquote(urlparse(dest).path), re.I)
            if m is None:
                continue
            ident = m[1]
            obj = found.setdefault(ident, {"id": ident, "aliases": [], "position": position, "pages": []})
            obj["aliases"] = dedup(obj["aliases"] + [label, ident])
            obj["pages"] = dedup(obj["pages"] + [dest])
            if separation(obj["position"], position) > 0.03:
                obj["conflicting_positions"] = True
    for obj in found.values():
        obj["pages"].sort(key=lambda u: (0 if "/star_data/" in u else 1 if "/v7.7l/" in u else 2, u))
        obj["page"] = obj["pages"][0]
    return list(found.values()), []


def separation(a, b):
    dra, ddec = (math.radians(x - y) for x, y in zip(a, b))
    hav = math.sin(ddec/2)**2 + math.cos(math.radians(a[1])) * math.cos(math.radians(b[1])) * math.sin(dra/2)**2
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, max(0.0, hav)))))


def coarse_position(ident):
    m = re.match(r"J(\d\d)(\d\d)([+\-])(\d\d)(\d)", ident)
    if not m:
        return None
    return [15 * (int(m[1]) + (int(m[2]) + 0.5)/60), (1 if m[3] == "+" else -1) *
            (int(m[4]) + (int(m[5]) + 0.5)/10)]


def crossmatch(rows, maxi):
    alias_map = {}
    for m in maxi:
        for alias in m["aliases"]:
            if key(alias):
                alias_map.setdefault(key(alias), set()).add(m["id"])
    by_id = {m["id"]: m for m in maxi}
    inverse = {}
    for b in rows:
        aliases = set().union(*(alias_map.get(key(a), set()) for a in b["aliases"]))
        coords = set()
        if b["position"]:
            coords = {m["id"] for m in maxi if m["position"] and separation(m["position"], b["position"]) <= 1/60}
        choices = sorted(aliases | coords)
        b.update(maxi_id=None, maxi_name=None, match="unmatched", match_candidates=choices)
        if len(choices) > 1:
            b["match"] = "ambiguous"
            continue
        if not choices:
            continue
        ident = choices[0]
        m = by_id[ident]
        # A MAXI J-name is a coarse sanity check only. It never supplies a
        # positional match. Measured catalogue coordinates can supply one.
        p = m["position"] or coarse_position(ident)
        if m.get("conflicting_positions"):
            b["match"] = "maxi_conflicting_positions"
            continue
        if p and b["position"] and separation(p, b["position"]) > (0.5 if not m["position"] else 0.05):
            b["match"] = "position_conflict"
            continue
        b.update(maxi_id=ident, maxi_name=next((a for a in m["aliases"] if a != ident), ident),
                 maxi_page=m["page"], maxi_pages=m["pages"], match="exact_alias" if ident in aliases else "position_within_1arcmin")
        inverse.setdefault(ident, []).append(b)
    for group in inverse.values():
        if len(group) > 1:
            for b in group:
                b.update(maxi_id=None, maxi_name=None, match="nonunique_reverse")


def official_roster(initial, parser, label):
    pending = list(initial)
    visited, rows, provenance = set(), {}, []
    while pending:
        url = pending.pop(0)
        if url.rstrip("/") in visited:
            continue
        visited.add(url.rstrip("/"))
        if len(visited) > 40:
            raise ValueError(f"unexpected {label} navigation; refusing partial inventory")
        try:
            data, meta = fetch(url)
        except Exception as exc:
            provenance.append({"url": url, "error": str(exc)})
            continue
        new, more = parser(data, meta["url"])
        meta["source_links"] = len(new)
        provenance.append(meta)
        for obj in new:
            prior = rows.get(obj["id"])
            if prior:
                prior["aliases"] = dedup(prior["aliases"] + obj["aliases"])
                prior["position"] = prior["position"] or obj["position"]
            else:
                rows[obj["id"]] = obj
        pending.extend(more)
    return sorted(rows.values(), key=lambda x: x["id"]), provenance


def save_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=1, sort_keys=True) + "\n")


def inventory(out):
    bat, bat_log = official_roster([BAT], bat_index, "BAT")
    maxi, maxi_log = official_roster([urljoin(MAXI, "top/slist.html")], maxi_index, "MAXI")
    crossmatch(bat, maxi)
    for name, url in (("maxi-lc-readme.txt", urljoin(MAXI, "top/lc_readme.txt")),):
        try:
            doc, meta = fetch(url)
            out.mkdir(parents=True, exist_ok=True)
            (out/name).write_bytes(doc)
            maxi_log.append(dict(meta, description="archive column documentation"))
        except requests.RequestException as exc:
            maxi_log.append({"url":url, "error":str(exc)})
    save_json(out/"inventory.json", {"sources": bat, "maxi_sources": maxi, "bat_indexes": bat_log,
                                      "maxi_indexes": maxi_log, "inventory_unix": int(time.time())})
    print(f"BAT roster {len(bat)}, MAXI roster {len(maxi)}, matched {sum(bool(b['maxi_id']) for b in bat)}", flush=True)
    if len(bat) < 900 or not any("/weak/" in b["bat_page"] for b in bat):
        raise RuntimeError("BAT roster coverage check failed (see inventory.json)")
    if len(maxi) < 150 or sum(bool(b["maxi_id"]) for b in bat) < 75:
        raise RuntimeError("MAXI roster or crossmatch coverage check failed (see inventory.json)")


def normalized_column(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def ascii_curve(data, instrument):
    txt = data.decode("utf-8", errors="replace")
    if "<html" in txt[:1000].lower() or "<!doctype" in txt[:1000].lower():
        raise ValueError("received HTML instead of light curve")
    comments, lines = [], []
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        if line[0] in "#!;" or not re.match(r"^[+\-.\d]", line):
            comments.append(line)
        else:
            lines.append(line)
    if not lines:
        raise ValueError("no numeric data rows")
    a = np.loadtxt(io.StringIO("\n".join(lines)), ndmin=2)
    if a.shape[1] < 3:
        raise ValueError("expected MJD, rate and uncertainty")
    flag = None
    if instrument == "bat":
        # Daily BAT ASCII includes DATA_FLAG, if present, as a named scalar
        # column. Reject nonzero flag values. FITS is preferred.
        for comment in comments:
            fields = [normalized_column(s) for s in comment.lstrip("#!; ").split()]
            if len(fields) == a.shape[1] and "dataflag" in fields:
                flag = a[:, fields.index("dataflag")]
                break
    return a[:, 0], a[:, 1], a[:, 2], flag


def bat_fits(data):
    with fits.open(io.BytesIO(data), memmap=False) as hdus:
        hdu = next(h for h in hdus if isinstance(h, fits.BinTableHDU) and h.data is not None)
        names = {normalized_column(c): c for c in hdu.columns.names}
        def col(*options):
            return next(names[c] for c in options if c in names)
        tname, rname, ename = col("time", "mjd"), col("rate", "flux"), col("error", "rateerr", "rateerror", "err")
        table = hdu.data
        t = np.array(table[tname], dtype=float)
        h = hdu.header
        unit = (hdu.columns[tname].unit or h.get("TIMEUNIT") or "d").lower()
        if unit in {"s", "sec", "second", "seconds"}:
            zero = float(h.get("MJDREF", h.get("MJDREFI", 0) + h.get("MJDREFF", 0)))
            if not zero:
                raise ValueError("time in seconds missing MJDREF")
            t = zero + (t + float(h.get("TIMEZERO", 0))) / 86400.
        elif unit not in {"d", "day", "days", "mjd"}:
            raise ValueError(f"unrecognized BAT FITS time unit {unit}")
        elif np.nanmedian(t) < 40000:
            t = t + float(h.get("MJDREF", h.get("MJDREFI", 0) + h.get("MJDREFF", 0))) + float(h.get("TIMEZERO", 0))
        flagname = names.get("dataflag") or names.get("quality")
        flag = np.array(table[flagname]) if flagname else None
        return t, np.array(table[rname], dtype=float), np.array(table[ename], dtype=float), flag


def clean_curve(raw):
    t, y, dy, flag = raw
    if not (t.ndim == y.ndim == dy.ndim == 1):
        raise ValueError("expected single-band daily series")
    good = np.isfinite(t) & np.isfinite(y) & np.isfinite(dy) & (dy > 0) & (t > 50000) & (t < 85000)
    if flag is not None:
        good &= (flag == 0)
    t, y, dy = (q[good] for q in (t, y, dy))
    if not len(t):
        raise ValueError("no finite daily observations with positive errors and good flags")
    idx = np.argsort(t)
    t, y, dy = t[idx], y[idx], dy[idx]
    if np.any(np.diff(t) == 0):
        raise ValueError("duplicate time stamps in daily archive product")
    return t, y, dy


def candidates(row, instrument):
    if instrument == "bat":
        page = row["bat_page"]
        base = page.rstrip("/")
        initial = [base + suffix for suffix in (".lc.fits", ".lc.txt")]
        pages = [page]
    else:
        ident = row["maxi_id"]
        pages = row["maxi_pages"]
        initial = [urljoin(p, ident + "_g_lc_1day_all.dat") for p in pages if "/star_data/" in p]
    yield from initial
    for page in pages:
        try:
            html, meta = fetch(page)
        except requests.RequestException:
            continue
        ranked = []
        for _, dest, label in page_links(html, meta["url"]):
            path = unquote(urlparse(dest).path.lower())
            if instrument == "bat":
                if "orbit" in path or not re.search(r"\.(fits|txt)(\.gz)?$", path):
                    continue
                if not ("daily" in label.lower() or ".lc." in path or "day" in path):
                    continue
                rank = 0 if ".fits" in path else 1
            else:
                if "ssc" in path or "orb" in path or not re.search(r"\.(dat|txt)(\.gz)?$", path):
                    continue
                if not re.search(r"1day|1d_|day_", path) or not ("all" in path or "all" in label.lower()):
                    continue
                rank = 0 if "gsc" in path or "_g_" in path else 1
            ranked.append((rank, dest))
        for _, dest in sorted(ranked):
            if dest not in initial:
                yield dest


def load_curve(row, instrument, out):
    result, attempts = {"instrument": instrument, "status": "unavailable"}, []
    result["attempts"] = attempts
    if instrument == "maxi" and not row.get("maxi_id"):
        return {"instrument": instrument, "status": "no_unique_counterpart", "reason": row["match"]}
    try:
        for url in dict_order(candidates(row, instrument)):
            try:
                data, meta = fetch(url)
                if data.startswith(b"\x1f\x8b"):
                    data = gzip.decompress(data)
                raw = bat_fits(data) if instrument == "bat" and data.startswith(b"SIMPLE") else ascii_curve(data, instrument)
                t, y, dy = clean_curve(raw)
                p = out/"curves"/f"{row['id']}.{instrument}.npz"
                np.savez_compressed(p, t=t, y=y, dy=dy)
                result.update(meta, status="ok", rows_raw=len(raw[0]), rows_used=len(t),
                              first_mjd=float(t[0]), last_mjd=float(t[-1]),
                              negative_rates=int(sum(y < 0)), quality_column=raw[3] is not None,
                              band="15-50 keV" if instrument == "bat" else "2-20 keV")
                return result
            except Exception as exc:
                attempts.append({"url": url, "error": str(exc)[:320]})
    except Exception as exc:
        attempts.append({"error": str(exc)[:320]})
    return result


def dict_order(it):
    """Deduplicate lazily so successful direct URLs never need an HTML fetch."""
    seen = set()
    for item in it:
        if item not in seen:
            seen.add(item)
            yield item


def curves(infile, out, shard, shards, workers):
    if not 0 <= shard < shards or workers < 1:
        raise ValueError("invalid shard or worker count")
    roster = json.loads(infile.read_text())["sources"][shard::shards]
    (out/"curves").mkdir(parents=True, exist_ok=True)
    def one(row):
        return {"id": row["id"], "name": row["name"],
                **{k: load_curve(row, k, out) for k in ("bat", "maxi")}}
    with (out/f"ingest-{shard}.jsonl").open("w") as trace, concurrent.futures.ThreadPoolExecutor(workers) as pool:
        for rec in pool.map(one, roster):
            trace.write(json.dumps(rec, sort_keys=True) + "\n")
            trace.flush()
    print(f"finished {len(roster)} inventory rows in shard {shard}/{shards}", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="mode", required=True)
    ip = subs.add_parser("inventory")
    ip.add_argument("--out", type=Path, required=True)
    cp = subs.add_parser("curves")
    cp.add_argument("--inventory", type=Path, required=True)
    cp.add_argument("--out", type=Path, required=True)
    cp.add_argument("--shard", type=int, default=0)
    cp.add_argument("--shards", type=int, default=1)
    cp.add_argument("--workers", type=int, default=3)
    args = p.parse_args()
    if args.mode == "inventory":
        inventory(args.out)
    else:
        curves(args.inventory, args.out, args.shard, args.shards, args.workers)
