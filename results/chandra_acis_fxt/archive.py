"""Inventory released ACIS pointings and acquire archive-native Level 2 bytes.

CXC schema and archive layout
https://cda.harvard.edu/cxctap/tables
https://cxc.cfa.harvard.edu/ciao/ahelp/download_chandra_obsid.html
https://cxc.cfa.harvard.edu/ciao/download/sciserver.html

The independent census is by public_avail_date. Observation start time, target
class, pointing, exposure length and measured source counts are never cuts.
"""
import argparse
import datetime as dt
import gzip
import hashlib
import http.client
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path

TAP = "https://cda.cfa.harvard.edu/cxctap/sync"
ROOTS = (
    "https://cxc.cfa.harvard.edu/cdaftp/byobsid",
    "https://heasarc.gsfc.nasa.gov/FTP/chandra/data/byobsid",
)
AGENT = "Open-Space-Chandra-FXT/1 (github.com/dovvnloading/Open-Space)"
CHUNK = 1024 * 1024
FIELDS = (
    "obsid", "instrument", "grating", "exposure_mode", "status",
    "public_avail_date", "start_date", "exposure_time", "event_count",
    "target_name", "ra", "dec", "proposal_number",
)
PROPOSAL_FIELDS = ("proposal_type", "proposal_title", "proposal_category", "proposal_pi")
ARCHIVE_HOSTS = {
    "cda.harvard.edu", "cda.cfa.harvard.edu", "cxc.harvard.edu",
    "cxc.cfa.harvard.edu", "heasarc.gsfc.nasa.gov",
}


def dump(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def check_archive_url(url):
    location = urllib.parse.urlsplit(url)
    if (location.scheme != "https" or location.hostname not in ARCHIVE_HOSTS
            or location.username is not None or location.password is not None
            or location.port not in {443, None}):
        raise ValueError(f"unexpected archive URL: {url}")


class ArchiveRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_archive_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download(url, dest, data=None):
    """Write the HTTP entity verbatim, validating transport length when provided."""
    check_archive_url(url)
    opener = urllib.request.build_opener(ArchiveRedirect())
    request = urllib.request.Request(
        url, data=data, headers={"User-Agent": AGENT, "Accept-Encoding": "identity"}
    )
    part = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(3):
        digest = hashlib.sha256()
        size = 0
        try:
            with opener.open(request, timeout=90) as response, part.open("wb") as out:
                while block := response.read(CHUNK):
                    out.write(block)
                    digest.update(block)
                    size += len(block)
                expected = response.headers.get("Content-Length")
                if size == 0 or expected is not None and int(expected) != size:
                    raise OSError(f"incomplete response: expected {expected}, got {size}")
                metadata = {
                    "url": url, "resolved_url": response.geturl(), "bytes": size,
                    "sha256": digest.hexdigest(), "last_modified": response.headers.get("Last-Modified"),
                    "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                }
            part.replace(dest)
            return metadata
        except (OSError, http.client.HTTPException) as error:
            part.unlink(missing_ok=True)
            if isinstance(error, urllib.error.HTTPError) and error.code in {400, 401, 403, 404, 410}:
                raise
            if attempt == 2:
                raise
            time.sleep(2 ** (attempt + 1))


def tap(outdir, label, query):
    """Fetch CXC TABLEDATA; retain exact raw XML compressed losslessly."""
    body = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "QUERY": query,
    }).encode("ascii")
    xml_path = outdir / (label + ".xml")
    source = download(TAP, xml_path, data=body)
    raw = xml_path.read_bytes()
    packed = outdir / (label + ".xml.gz")
    packed.write_bytes(gzip.compress(raw, mtime=0))
    xml_path.unlink()
    source.update({"query": query, "post_body": body.decode(), "raw_gzip": packed.name})
    root = ET.fromstring(raw)
    statuses = [i for i in root.findall(".//{*}INFO") if i.get("name", "").upper() == "QUERY_STATUS"]
    for status in statuses:
        if status.get("value", "").upper() not in {"OK", "OVERFLOW"}:
            raise ValueError(f"TAP {label}: {status.text or status.attrib}")
    table = root.find(".//{*}TABLE")
    if table is None or table.find(".//{*}TABLEDATA") is None:
        raise ValueError(f"TAP {label}: expected VOTable TABLEDATA")
    fields = [f.attrib["name"].lower().split(".")[-1] for f in table.findall("{*}FIELD")]
    rows = [dict(zip(fields, ((v.text or "").strip() for v in r), strict=True))
            for r in table.findall(".//{*}TR")]
    return rows, source


def inventory(outdir, start, end):
    outdir.mkdir(parents=True, exist_ok=True)
    first, last = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    if last < first:
        raise ValueError("end precedes start")
    # Paging by the integer primary key keeps coverage independent of TAP MAXREC
    # and avoids differences between TAP implementations' timestamp literals.
    predicate = "o.instrument LIKE 'ACIS%'"
    count_rows, count_source = tap(outdir, "census_count", (
        "SELECT COUNT(*) AS n FROM cxc.observation o WHERE " + predicate
    ))
    expected = int(count_rows[0]["n"])
    if expected < 1:
        raise ValueError("empty ACIS observation catalog")
    sources = [count_source]
    accepted = []
    exceptions = []
    inspected = page = 0
    lastid = -1
    while True:
        columns = ", ".join("o." + field for field in FIELDS)
        proposal_cols = ('p.type AS proposal_type, p.title AS proposal_title, '
                         'p.category AS proposal_category, p."PI" AS proposal_pi')
        query = (f"SELECT TOP 2000 {columns}, {proposal_cols} FROM cxc.observation o "
                 "LEFT OUTER JOIN cxc.proposal p ON p.number = o.proposal_number "
                 f"WHERE {predicate} AND o.obsid > {lastid} ORDER BY o.obsid")
        rows, source = tap(outdir, f"census_{page:03d}", query)
        sources.append(source)
        page += 1
        if not rows:
            break
        for row in rows:
            if set(row) != set(FIELDS + PROPOSAL_FIELDS):
                raise ValueError(f"unexpected CXC fields: {list(row)}")
            obsid = int(row["obsid"])
            if obsid <= lastid:
                raise ValueError(f"duplicate or out-of-order CXC ObsID {obsid}")
            lastid = obsid
            inspected += 1
            status = row["status"].upper()
            mode = row["exposure_mode"].upper()
            grating = row["grating"].upper()
            if grating in {"HETG", "LETG"} or mode in {"CC", "CONTINUOUS", "CONTINUOUS CLOCKING"}:
                continue
            release = row["public_avail_date"]
            if not release:
                # Retain catalog gaps affecting the interval. Observations predating
                # the interval can become public during it, even years later.
                if status == "ARCHIVED":
                    exceptions.append({**row, "reason": "missing_public_date"})
                continue
            date = dt.date.fromisoformat(release[:10])
            if not first <= date <= last:
                continue
            if status != "ARCHIVED":
                exceptions.append({**row, "reason": "not_marked_archived"})
                continue
            if grating not in {"NONE", "", "NO GRATING"} or mode not in {
                "TE", "TIMED", "TIMED EXPOSURE", "",
            }:
                # Do not discard unexpected labels. Confirm actual science
                # configuration from READMODE/GRATING of the downloaded EVENTS.
                exceptions.append({**row, "reason": "check_mode_from_FITS"})
            proposal_number = row["proposal_number"]
            proposal_link = None
            if proposal_number:
                escaped = proposal_number.replace("'", "''")
                query_text = ("SELECT * FROM cxc.proposal p WHERE p.number = "
                              f"'{escaped}'")
                proposal_link = TAP + "?" + urllib.parse.urlencode({
                    "REQUEST": "doQuery", "FORMAT": "text", "LANG": "ADQL", "QUERY": query_text,
                })
            accepted.append({
                "target_name": row["target_name"], "proposal_number": proposal_number,
                **{field: row[field] for field in PROPOSAL_FIELDS},
                "catalog_ra_deg": float(row["ra"]) if row["ra"] else None,
                "catalog_dec_deg": float(row["dec"]) if row["dec"] else None,
                "observation_url": f"https://cda.harvard.edu/chaser/startViewer.do?menuItem=details&obsid={obsid}",
                "proposal_url": proposal_link,
                "obsid": obsid, "instrument": row["instrument"], "grating": grating,
                "release_date": release, "start_date": row["start_date"],
                "exposure_mode": mode, "status": status,
                "exposure_s": 1000 * float(row["exposure_time"]) if row["exposure_time"] else None,
                "event_count": int(row["event_count"]) if row["event_count"] else None,
            })
    if inspected != expected:
        raise ValueError(f"CXC census changed or paginated incompletely: {inspected} of {expected}")
    if not accepted:
        raise ValueError("no archived ACIS imaging entries in interval; inspect census before interpreting")
    dump(outdir / "inventory_sources.json", {
        "sources": sources, "catalog_rows": expected, "selected_obsids": len(accepted),
        "start_inclusive": start, "end_inclusive": end,
        "source_sha256_definition": "Exact HTTP response entity, before local gzip compression",
    })
    dump(outdir / "inventory_exceptions.json", exceptions)
    dump(outdir / "inventory.json", accepted)
    return accepted


class Hrefs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = set()

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.links.update(v for k, v in attrs if k == "href" and v)


def product_names(raw_html, obsid):
    parser = Hrefs()
    parser.feed(raw_html)
    available = {}
    for href in parser.links:
        basename = urllib.parse.unquote(urllib.parse.urlsplit(href).path.rsplit("/", 1)[-1])
        match = re.fullmatch(r"acisf(\d+)(?:_[a-zA-Z0-9]+)*?N(\d+)([a-zA-Z0-9_]*)_evt2\.fits(?:\.gz)?", basename, re.I)
        if not match or int(match[1]) != obsid:
            continue
        key = re.sub(r"N\d+", "N#", basename.removesuffix(".gz"), count=1, flags=re.I)
        # Separate alternating-exposure streams. Ignore older revisions and
        # duplicate compressed/uncompressed representations of the same file.
        rank = (int(match[2]), basename.endswith(".gz"))
        if key not in available or rank > available[key][0]:
            available[key] = (rank, basename)
    return sorted(value[1] for value in available.values())


def inspect_product(path, obsid):
    from astropy.io import fits

    # Read the payload to EOF so gzip validates its CRC and length even when
    # the FITS reader only needs the EVENTS and GTI extensions.
    opener = gzip.open if path.name.endswith(".gz") else open
    payload = hashlib.sha256()
    with opener(path, "rb") as stream:
        while block := stream.read(CHUNK):
            payload.update(block)
    with fits.open(path, memmap=False) as hdus:
        event = hdus["EVENTS"]
        header = event.header
        if int(header["OBS_ID"]) != obsid or header["INSTRUME"].strip().upper() != "ACIS":
            raise ValueError("event product OBS_ID or INSTRUME mismatch")
        readmode = str(header.get("READMODE", "")).strip().upper()
        grating = str(header.get("GRATING", "")).strip().upper()
        if readmode != "TIMED" or grating != "NONE":
            raise ValueError(f"unverified imaging configuration GRATING={grating!r}, READMODE={readmode!r}")
        return {
            "fits_sha256": payload.hexdigest(), "rows": int(header["NAXIS2"]),
            "revision": header.get("REVISION"), "caldbver": header.get("CALDBVER"),
            "ascdsver": header.get("ASCDSVER"), "cycle": header.get("CYCLE"),
            "tstart": header.get("TSTART"), "tstop": header.get("TSTOP"),
            "readmode": readmode, "grating": grating,
            "ra_nom_deg": header.get("RA_NOM"), "dec_nom_deg": header.get("DEC_NOM"),
            "ra_targ_deg": header.get("RA_TARG"), "dec_targ_deg": header.get("DEC_TARG"),
        }


def fetch_events(obsid, dest):
    """Download all current, distinct timed-mode EVT2 streams for one ObsID."""
    obsid, dest = int(obsid), Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    errors = []
    attempts = []
    for idx, root in enumerate(ROOTS):
        base = f"{root}/{obsid % 10}/{obsid}/primary/"
        listing = dest / f"primary_{idx}.html"
        attempt = {"index_url": base, "downloaded": []}
        attempts.append(attempt)
        try:
            listing_source = download(base, listing)
            attempt["index"] = listing_source
            names = product_names(listing.read_text(errors="replace"), obsid)
            attempt["selected_files"] = names
            if not names:
                raise ValueError("no ACIS evt2 files in primary index")
            products = []
            for name in names:
                path = dest / name
                source = download(urllib.parse.urljoin(base, name), path)
                attempt["downloaded"].append(dict(source, filename=name))
                products.append({
                    **source, **inspect_product(path, obsid), "obsid": obsid, "filename": name,
                    "path": str(path), "index_url": base, "index_sha256": listing_source["sha256"],
                })
            dump(dest / "acquisition.json", attempts)
            dump(dest / "products.json", products)
            return products
        except (OSError, ValueError, KeyError, EOFError, http.client.HTTPException) as error:
            errors.append(f"{base}: {type(error).__name__}: {error}")
            attempt["error"] = errors[-1]
            dump(dest / "acquisition.json", attempts)
    raise RuntimeError(" | ".join(errors))


def replay(manifest, dest):
    """Retrieve exact archived products or fail explicitly on a revision change."""
    dest.mkdir(parents=True, exist_ok=True)
    for product in json.loads(manifest.read_text()):
        name = product["filename"]
        if Path(name).name != name:
            raise ValueError("unexpected path in product manifest")
        output = dest / name
        checked = dest / (name + ".verified")
        try:
            source = download(product["url"], checked)
            if (source["sha256"], source["bytes"]) != (product["sha256"], product["bytes"]):
                raise ValueError(f"{name} differs from manifest. Request original revision from CXC.")
            checked.replace(output)
        finally:
            checked.unlink(missing_ok=True)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    commands = cli.add_subparsers(dest="command", required=True)
    census = commands.add_parser("inventory")
    census.add_argument("outdir", type=Path)
    census.add_argument("--start", default="2021-01-01")
    census.add_argument("--end", default="2026-11-21")
    acquire = commands.add_parser("fetch")
    acquire.add_argument("obsid", type=int)
    acquire.add_argument("outdir", type=Path)
    pinned = commands.add_parser("replay")
    pinned.add_argument("manifest", type=Path)
    pinned.add_argument("outdir", type=Path)
    args = cli.parse_args()
    if args.command == "inventory":
        inventory(args.outdir, args.start, args.end)
    elif args.command == "fetch":
        fetch_events(args.obsid, args.outdir)
    else:
        replay(args.manifest, args.outdir)
