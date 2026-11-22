"""Positional counterpart checks. Preserve responses and replay them by hash.

Call vet_trigger() on real and scrambled triggers with the same criteria.
The output directory may be shared across triggers to reuse static inputs.

catalog_incomplete means a mandatory catalog failed response integrity,
schema, truncation or transport checks. known_transient records variability
previously published at the position in the same ObsID. foreground_veto
records a precise stellar/Galactic match. ambiguous includes uncertain
astrometry, competing matches, prior variability in other visits and
|b|<10 degrees. no_counterpart leaves the origin undetermined.

host_candidate requires a catalog galaxy with z>0.001 within max(5 arcsec,
the X-ray match radius), a one-sided upper chance-projection probability
<=0.01, and no unresolved competitors. That probability assumes a locally
homogeneous Poisson population at the catalog's depth. The 95% upper galaxy
density comes from the 30-120 arcsec annulus. It is not a posterior probability
of association. Catalogs can miss faint stars, nearby extended host disks
and off-nuclear events. A candidate still requires a literature/alert/pointing
and repeated-visit review to establish novelty.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import tempfile
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table
from scipy.stats import chi2

GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
SIMBAD = "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync"
NED = "https://ned.ipac.caltech.edu/tap/sync"
VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/votable"
MIRRORS = {
    "gaia": ("gaia_aip", "https://gaia.aip.de/tap/sync"),
    "simbad": ("simbad_idia", "https://simbad.idia.ac.za/simbad/sim-tap/sync"),
    "ned": ("ned_vo", "https://vo.ned.ipac.caltech.edu/tap/sync"),
    "quirola2023": ("quirola2023_idia", "https://vizier.idia.ac.za/viz-bin/votable"),
    "dillmann2025": ("dillmann2025_idia", "https://vizier.idia.ac.za/viz-bin/votable"),
}
CUDIDI = ("https://zenodo.org/records/19160904/files/"
          "final_candidates_catalogue_2026_clean_with_suffixes.csv?download=1")
CUDIDI_MD5 = "81e1fe86d8bb346ce1ca806f144a8d63"
CAMPAIGN_REF = "refs/heads/research/chandra-acis-fast-transients"
MAX_ROWS = 25000
FIELD = 120.0  # Arcseconds. Density annulus extends from FIELD/4 to FIELD.


def digest(data):
    return hashlib.sha256(data).hexdigest()


def save_bytes(path, data):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as staging:
        staging.write(data)
        temporary = Path(staging.name)
    temporary.replace(path)


def read_input(label, url, root, replay, records, params=None):
    if params:
        url += "?" + urlencode(params)
    key = digest(url.encode())[:24]
    name = f"{label}-{key}"
    manifest = root / f"{name}.json"
    if manifest.exists():
        memo = json.loads(manifest.read_text())
        if memo["url"] != url:
            raise ValueError("request cache key collision")
        if memo.get("ok"):
            body = (root / memo["body_file"]).read_bytes()
            if digest(body) != memo["sha256"]:
                raise ValueError(f"altered catalog response {memo['body_file']}")
            records[label] = memo
            return body
        if replay or memo.get("permanent_error"):
            records[label] = memo
            raise ValueError(f"{label} previously failed")
    if replay:
        raise ValueError(f"missing replay input {name}")
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != CAMPAIGN_REF
            or os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("Catalog retrieval requires the research branch's Actions campaign")
    memo = {"url": url, "attempts": [], "ok": False}
    for attempt in range(3):
        body = None
        status = None
        entry = {"time_utc": datetime.now(timezone.utc).isoformat()}
        try:
            request = Request(url, headers={"User-Agent": "Open-Space-ACIS-FXT/1", "Accept-Encoding": "identity"})
            with urlopen(request, timeout=65) as response:
                entry["final_url"] = response.geturl()
                status = response.status
                body = response.read()
        except HTTPError as error:
            status = error.code
            body = error.read()
            entry["error"] = str(error)
        except (URLError, TimeoutError, OSError) as error:
            entry["error"] = str(error)
        entry["http_status"] = status
        if body is not None:
            sha = digest(body)
            path = f"{name}-{sha[:20]}.bin"
            save_bytes(root / path, body)
            entry.update(body_file=path, sha256=sha, size_bytes=len(body))
        memo["attempts"].append(entry)
        if status == 200 and body:
            memo.update(ok=True, body_file=path, sha256=sha, size_bytes=len(body))
        elif status and 400 <= status < 500 and status not in {408, 425, 429}:
            memo["permanent_error"] = True
        save_bytes(manifest, (json.dumps(memo, indent=2) + "\n").encode())
        records[label] = memo
        if memo["ok"]:
            return body
        if memo.get("permanent_error"):
            break
        if attempt < 2:
            time.sleep(2 ** (attempt + 1))
    raise ValueError(f"{label} retrieval failed")


def catalog_rows(label, endpoint, query, required, root, replay, records, receipt_label=None):
    # Small published lists are cached whole and checked against the CDS ReadMe.
    fixed = {
        "quirola2023": ("J/A+A/675/A44/table3", "**", 8),
        "dillmann2025": ("J/MNRAS/537/931/catalog", "CSCName,RACdeg,DECdeg,CSCobsId,TransType", 3559),
    }.get(label)
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "MAXREC": MAX_ROWS, "QUERY": query}
    if fixed:
        params = {"-source": fixed[0], "-out": fixed[1], "-out.max": MAX_ROWS}
    if label == "quirola2023":
        params.update({"-out.add": "_RAJ2000,_DEJ2000", "-oc.form": "d"})
    receipt_label = receipt_label or label
    raw = read_input(receipt_label, endpoint, root, replay, records, params)
    decode = decode_published_table if fixed else decode_table
    rows = decode(raw, tuple(required), fixed[2] if fixed else None, label)
    records[receipt_label]["row_count"] = len(rows)
    return rows


def catalog_angle(value, unit, is_ra):
    text = str(value).strip()
    if unit in {"deg", "degree", "degrees"}:
        v = float(text)
    elif unit in ({"h:m:s", "hms"} if is_ra else {"d:m:s", "dms"}):
        parts = re.split(r"[:\s]+", text)
        if len(parts) != 3:
            raise ValueError("Incomplete sexagesimal coordinate")
        h, m, s = map(float, parts)
        if not (0 <= m < 60 and 0 <= s < 60):
            raise ValueError("Invalid sexagesimal coordinate")
        v = (abs(h) + m / 60 + s / 3600) * (-1 if text.startswith("-") else 1) * (15 if is_ra else 1)
    else:
        raise ValueError(f"Unexpected coordinate unit {unit!r}")
    if not (0 <= v < 360 if is_ra else -90 <= v <= 90):
        raise ValueError("Invalid catalog coordinate")
    return v


def decode_table(raw, required, expected_count, label=None):
    doc = ET.fromstring(raw)
    states = [el for el in doc.iter() if el.tag.rsplit("}", 1)[-1] == "INFO"
              and el.get("name", "").upper() == "QUERY_STATUS"]
    if (not states and expected_count is None) or any(el.get("value", "").upper() != "OK" for el in states):
        raise ValueError("TAP status missing/overflow/error " + " ".join(str(el.text)[:200] for el in states))
    vo = parse_single_table(io.BytesIO(raw))
    table = vo.to_table(use_names_over_ids=True)
    columns = {str(n).lower(): n for n in table.colnames}
    fields = {str(field.name).lower(): field for field in vo.fields}
    if label == "ned":
        types = [k for k in ("preftype", "prefobjtype", "pref_objtype", "pref_obj_type", "preferred_type", "pretype",
                              "objtype", "obj_type", "object_type", "type") if k in columns]
        if not types:
            types = [k for k, f in fields.items() if (f.ucd or "").lower().split(";")[0] == "src.class"
                     or "preferred object type" in (f.description or "").lower()]
        if len(types) != 1:
            raise ValueError(f"NED object-type field is missing or ambiguous: {sorted(columns)}")
        columns["pretype"] = columns[types[0]]
    if label == "quirola2023":
        pairs = [("_raj2000","_dej2000"), ("_ra.icrs","_de.icrs"), ("_ra.icrs","_dec.icrs"),
                 ("radeg","dedeg"), ("raj2000","dej2000"), ("raicrs","deicrs")]
        pair = next((p for p in pairs if all(k in columns for k in p)), None)
        if pair is None:
            raise ValueError(f"VizieR coordinate pair missing: {sorted(columns)}")
        columns.update(radeg=columns[pair[0]], dedeg=columns[pair[1]])
        xml_units = {e.get("name", "").lower(): e.get("unit", "").strip().lower() for e in doc.iter()
                     if e.tag.rsplit("}", 1)[-1] == "FIELD"}
        units = {key: xml_units[key] for key in pair}
    if (not set(required) <= columns.keys() or len(table) >= MAX_ROWS or
            (expected_count is not None and len(table) != expected_count)):
        raise ValueError("Incomplete schema or unexpected catalog row count")
    def scalar(x):
        if getattr(x, "mask", False):
            return None
        if isinstance(x, bytes):
            return x.decode("utf8").strip()
        return str(x).strip()
    rows = [{k: scalar(row[columns[k]]) for k in required} for row in table]
    if label == "ned" and rows and not any(row["pretype"] for row in rows):
        raise ValueError("No reported object types in a nonempty NED cone")
    if label == "quirola2023":
        for row in rows:
            row["radeg"] = catalog_angle(row["radeg"], units[pair[0]], True)
            row["dedeg"] = catalog_angle(row["dedeg"], units[pair[1]], False)
    return rows


decode_published_table = lru_cache(maxsize=2)(decode_table)


def number(value):
    try:
        val = float(value)
        return val if math.isfinite(val) else None
    except (ValueError, TypeError):
        return None


def separation(ra, dec, ra2, dec2):
    # Numerically stable at small separations, across the meridian and at the poles.
    if not (0 <= ra2 <= 360 and -90 <= dec2 <= 90):
        raise ValueError("catalog position missing or out of range")
    d1, d2 = math.radians(dec), math.radians(dec2)
    h = (math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) *
         math.sin(math.radians(ra2 - ra) / 2) ** 2)
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, max(0.0, h))))) * 3600


def pm_position(row, epoch):
    """First-order great-circle propagation. Gaia pmra includes cos(declination)."""
    ra, dec = float(row["ra"]), float(row["dec"])
    pmra, pmdec = number(row["pmra"]), number(row["pmdec"])
    if pmra is None or pmdec is None or epoch is None:
        return ra, dec
    dt = epoch - float(row["ref_epoch"])
    dx, dy = pmra * dt / 1000, pmdec * dt / 1000
    bearing, step = math.atan2(dx, dy), math.radians(math.hypot(dx, dy) / 3600)
    lat = math.radians(dec)
    along = math.cos(lat) * math.cos(step) - math.sin(lat) * math.sin(step) * math.cos(bearing)
    east = math.sin(step) * math.sin(bearing)
    north = math.sin(lat) * math.cos(step) + math.cos(lat) * math.sin(step) * math.cos(bearing)
    return (ra + math.degrees(math.atan2(east, along))) % 360, math.degrees(math.atan2(north, math.hypot(along, east)))


def chance_upper(sep, count):
    """95% upper Poisson surface-density estimate from the 30-120 arcsec annulus."""
    upper = 0.5 * chi2.ppf(0.95, 2 * count + 2)
    return float(-math.expm1(-upper * sep ** 2 / (FIELD ** 2 - (FIELD / 4) ** 2)))


def cone(ra, dec, radius, x="ra", y="dec", frame="ICRS"):
    return (f"CONTAINS(POINT('{frame}',{x},{y}),"
            f"CIRCLE('{frame}',{ra:.10f},{dec:.10f},{radius/3600:.10f}))=1")


def is_star(row):
    hierarchy = {s.strip() for s in (row.get("path") or "").split(">")}
    return ("*" in hierarchy or "**" in hierarchy or row.get("otype") in
            {"*", "**", "XB", "HXB", "LXB", "Psr", "CV*", "SNR", "PWN", "GlC", "OpC"})


def is_host(label, row):
    # Exclude identified active nuclei; their variability has a competing explanation.
    if label == "ned":
        return row["pretype"] == "G"
    return row["otype"] in {"G", "GiG", "GiP", "LSB", "EmG", "H2G", "SBG"}


def gaia(row, dt):
    def sn(k):
        err = number(row[k + "_error"])
        return (number(row[k]) or 0) / err if err and err > 0 else 0
    high = sn("parallax") >= 5 or max(abs(sn("pmra")), abs(sn("pmdec"))) >= 5
    ruwe = number(row["ruwe"])
    foreground = high and ruwe is not None and ruwe <= 1.4
    # Bound coordinate+proper-motion uncertainties even with correlated solutions.
    err = max(number(row.get("ra_error")) or 0, number(row.get("dec_error")) or 0)
    pmerr = max(number(row.get("pmra_error")) or 0, number(row.get("pmdec_error")) or 0)
    uncertainty = 3 * math.sqrt(2) * (err + abs(dt) * pmerr) / 1000
    # Propagated barycentric coordinates omit annual parallax. Permit its full amplitude.
    uncertainty += max(0, number(row.get("parallax")) or 0) / 1000
    return foreground, uncertainty


def vet_trigger(trigger, output_dir, *, replay=False):
    """Return screening evidence. A blank optical field does not identify a host."""
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    ra, dec = float(trigger["ra_deg"]), float(trigger["dec_deg"])
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError("invalid ICRS position")
    epoch = number(trigger.get("epoch_year"))
    if epoch is None and number(trigger.get("event_mjd")) is not None:
        epoch = 2000.0 + (float(trigger["event_mjd"]) - 51544.5) / 365.25
    dt = (epoch - 2016.0) if epoch is not None else 0
    r90 = number(trigger.get("position_error_arcsec"))
    localization_valid = r90 is not None and 0 < r90 <= 15
    radius = max(1.0, math.sqrt(2) * (r90 if r90 and r90 > 0 else 5))
    bg = float(SkyCoord(ra=ra*u.deg, dec=dec*u.deg).galactic.b.deg)
    pos = {"ra_deg": ra, "dec_deg": dec, "galactic_latitude_deg": bg,
           "epoch_year": epoch, "position_error90_arcsec": r90, "match_radius_arcsec": radius}
    records, errors, matches, host_options = {}, {}, [], []
    # Include a 12 arcsec/yr envelope for high-proper-motion sources before propagation.
    # Cache an enclosing tile cone. Matching and density estimates retain the
    # original centroid. A half-step is 15" in each coordinate, hence <22".
    grid_ra, grid_dec = (round(ra*120)/120) % 360, round(dec*120)/120
    grid_field = FIELD + 22
    gaia_radius = grid_field + 30 * math.ceil(max(radius, 12 * abs(dt)) / 30)
    columns = ("source_id,ra,dec,ref_epoch,pmra,pmdec,pmra_error,pmdec_error,"
               "ra_error,dec_error,parallax,parallax_error,ruwe,phot_g_mean_mag")
    specs = {
        "gaia": (GAIA, f"SELECT {columns} FROM gaiadr3.gaia_source WHERE {cone(grid_ra,grid_dec,gaia_radius)}",
                 columns.split(",")),
        "simbad": (SIMBAD, "SELECT b.main_id,b.ra,b.dec,b.otype,b.rvz_redshift,"
                   "b.coo_err_maj,t.path,t.is_candidate FROM basic AS b "
                   "LEFT OUTER JOIN otypedef AS t ON b.otype=t.otype WHERE " +
                   cone(grid_ra,grid_dec,grid_field,"b.ra","b.dec"),
                   ["main_id","ra","dec","otype","rvz_redshift","path","is_candidate","coo_err_maj"]),
        "ned": (NED, "SELECT * FROM objdir WHERE " +
                cone(grid_ra,grid_dec,grid_field,frame="J2000"), ["prefname","ra","dec","pretype","z","uncmaja"]),
        "quirola2023": (VIZIER, None, ["xrt","radeg","dedeg","obsid","epos"]),
        "dillmann2025": (VIZIER, None, ["cscname","racdeg","decdeg","cscobsid","transtype"]),
    }
    for label, (endpoint, query, required) in specs.items():
        try:
            try:
                rows = catalog_rows(label, endpoint, query, required, root, replay, records)
            except (OSError, ValueError, LookupError, TypeError, ET.ParseError) as exc:
                if label not in MIRRORS:
                    raise
                records.setdefault(label, {})["query_error"] = str(exc)[:400]
                receipt_label, mirror = MIRRORS[label]
                rows = catalog_rows(label, mirror, query, required, root, replay, records, receipt_label)
            nearby = []
            for row in rows:
                if label == "gaia":
                    x, y = pm_position(row, epoch)
                elif label == "quirola2023":
                    x, y = float(row["radeg"]), float(row["dedeg"])
                elif label == "dillmann2025":
                    x, y = float(row["racdeg"]), float(row["decdeg"])
                else:
                    x, y = float(row["ra"]), float(row["dec"])
                dist = separation(ra, dec, x, y)
                if dist <= FIELD:
                    nearby.append((row, dist, x, y))
            disk = [(r, d) for r, d, _, _ in nearby if FIELD/4 < d <= FIELD]
            n_hosts = sum(is_host(label, r) for r, _ in disk) if label in {"simbad","ned"} else 0
            for row, sep, x, y in nearby:
                limit, kind = radius, "unclassified"
                if label == "gaia":
                    foreground, err = gaia(row, dt)
                    limit = math.hypot(radius, err)
                    if foreground:
                        kind = "stellar_astrometry"
                elif label in {"quirola2023", "dillmann2025"}:
                    limit = math.hypot(radius, max(2, number(row.get("epos")) or 0))
                    kind = "published_variability"
                elif label == "simbad":
                    limit = math.hypot(radius, 3 * (number(row["coo_err_maj"]) or 0) / 1000)
                    if is_star(row):
                        kind = ("stellar_or_galactic" if str(row["is_candidate"]).lower() in {"0","false"}
                                else "possible_stellar_or_galactic")
                elif label == "ned":
                    limit = math.hypot(radius, 1.25 * (number(row["uncmaja"]) or 0))
                    if row["pretype"] in {"*", "**", "V*", "WD*", "Psr", "Nova"}:
                        kind = "stellar_or_galactic"
                host = label in {"simbad","ned"} and is_host(label, row)
                if sep > max(limit, 10 if host else 0):
                    continue
                name = next(row[k] for k in ("source_id","main_id","prefname","xrt","cscname") if k in row)
                item = {"catalog": label, "name": name, "ra_deg": x, "dec_deg": y,
                        "separation_arcsec": sep, "radius_arcsec": limit, "kind": kind,
                        "chance_upper95": chance_upper(max(limit,sep), n_hosts if host else len(disk)), "data": row}
                if host:
                    redshift = number(row.get("z") or row.get("rvz_redshift"))
                    item.update(kind="galaxy", redshift=redshift)
                    host_options.append(item)
                matches.append(item)
        except (OSError, ValueError, TypeError, LookupError, ET.ParseError) as exc:
            errors[label] = str(exc)[:400]

    try:
        data = read_input("cudidi2026", CUDIDI, root, replay, records)
        if hashlib.md5(data, usedforsecurity=False).hexdigest() != CUDIDI_MD5:
            raise ValueError("file differs from the publisher's version-3 checksum")
        reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
        if not {"OBSID","RA","DEC","CANDIDATE"} <= set(reader.fieldnames or []):
            raise ValueError("CuDiDi schema mismatch")
        for row in reader:
            x, y = number(row["RA"]), number(row["DEC"])
            if x is None or y is None:
                coord = SkyCoord(row["RA"], row["DEC"], unit=(u.hourangle,u.deg))
                x, y = float(coord.ra.deg), float(coord.dec.deg)
            sep = separation(ra, dec, x, y)
            limit = math.hypot(radius, 3)
            if sep <= limit:
                matches.append({"catalog": "cudidi2026", "name": row["CANDIDATE"],
                                "kind": "published_variability", "separation_arcsec": sep,
                                "radius_arcsec": limit, "ra_deg": x, "dec_deg": y,
                                "data": {"obsid": row["OBSID"], **row}})
    except (OSError, ValueError, KeyError, TypeError) as exc:
        errors["cudidi2026"] = str(exc)[:400]

    prior = [m for m in matches if m["kind"] == "published_variability"]
    same_obs = [m for m in prior if number(trigger.get("obsid")) is not None and
                number(m["data"].get("obsid") or m["data"].get("cscobsid")) == number(trigger["obsid"])
                and m["data"].get("transtype") != "D"]
    foreground = [m for m in matches if m["kind"] in {"stellar_astrometry","stellar_or_galactic"}
                  and m["chance_upper95"] <= 0.01]
    uncertain = [m for m in matches if m["kind"] not in {"galaxy","published_variability"}]
    host_options.sort(key=lambda m: m["separation_arcsec"])
    viable = [m for m in host_options if m["redshift"] is not None and m["redshift"] > 0.001
              and m["chance_upper95"] <= 0.01 and m["separation_arcsec"] <= max(radius, 5)]
    host = viable[0] if viable else None
    competing = bool(host and any(separation(host["ra_deg"],host["dec_deg"],m["ra_deg"],m["dec_deg"]) > 1
                                  for m in host_options))
    if errors:
        status, reason = "catalog_incomplete", "Required input missing, truncated or unreadable"
    elif same_obs:
        status, reason = "known_transient", "Published variability at this position in this ObsID"
    elif epoch is None or not localization_valid:
        status, reason = "ambiguous", "Observation epoch or usable positional uncertainty missing"
    elif foreground:
        status, reason = "foreground_veto", "Close astrometric/typed stellar or Galactic counterpart"
    elif prior or uncertain or competing or abs(bg) < 10:
        status, reason = "ambiguous", "Prior variability, unclassified counterpart, host confusion or low Galactic latitude"
    elif host:
        status, reason = "host_candidate", "Single candidate galaxy with catalog redshift and a low chance projection estimate"
    elif host_options:
        status, reason = "ambiguous", "Nearby galaxy association fails the adopted host criteria"
    else:
        status, reason = "no_counterpart", "No catalog counterpart in the screening region"
    return {"status": status, "reason": reason, "position": pos, "host": host,
            "matches": sorted(matches, key=lambda m: m["separation_arcsec"]), "queries": records,
            "errors": errors, "coverage": "Published Chandra event lists: Quirola-Vasquez 2023, "
            "Dillmann 2025, Saathoff 2026. Catalog absence alone cannot establish a new event or an "
            "extragalactic origin. Galaxy redshifts refer to the galaxy. Large or offset hosts can be missed."}


def contract_probe(output_dir):
    """Exercise real services before fan-out using a frozen, instrument-vetoed trigger."""
    run = os.environ["GITHUB_RUN_ID"] + "-" + os.environ["GITHUB_RUN_ATTEMPT"]
    root = Path(output_dir)
    cache = root / run
    cache.mkdir(parents=True, exist_ok=True)
    source = {"obsid": 31579, "trigger_key": "real/0", "ra_deg": 350.89949962074587,
              "dec_deg": 58.80101615192483, "event_mjd": 61352.65428720553,
              "position_error_arcsec": 1.3394525120377347,
              "input_sha256": "aa401202b7c0ce30c70880f60964fe1584bc5604a77931817fe05b0d0782f865"}
    report = {"run": run, "cache": run, "source": source, "status": "incomplete"}
    try:
        result = vet_trigger(source, cache)
        report.update(queries=result["queries"], errors=result["errors"])
        if result["errors"]:
            raise ValueError("Service validation failed: " + json.dumps(result["errors"]))
        required = {"gaia": None, "simbad": None, "ned": None, "quirola2023": 8, "dillmann2025": 3559}
        for name, count in required.items():
            routes = [name] + ([MIRRORS[name][0]] if name in MIRRORS else [])
            if not any(isinstance((r := result["queries"].get(route, {})).get("row_count"), int)
                       and r.get("ok") and r.get("sha256") and (count is None or r["row_count"] == count)
                       for route in routes):
                raise ValueError("Missing validated " + name + " receipt")
        csv = result["queries"].get("cudidi2026", {})
        if not (csv.get("ok") and csv.get("sha256")):
            raise ValueError("Missing validated CuDiDi receipt")
        report["status"] = "validated"
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        save_bytes(root / "contract.json", (json.dumps(report, indent=2) + "\n").encode())
    print("Validated Gaia, SIMBAD, NED and three published Chandra transient lists.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trigger_json", type=Path, nargs="?")
    parser.add_argument("output_dir", type=Path, nargs="?")
    parser.add_argument("--probe", type=Path, metavar="DIR", help="run the real-catalog contract and retain its inputs")
    parser.add_argument("--replay", action="store_true", help="require all saved inputs; no network")
    args = parser.parse_args()
    if args.probe is not None:
        contract_probe(args.probe)
    elif args.trigger_json and args.output_dir:
        print(json.dumps(vet_trigger(json.loads(args.trigger_json.read_text()), args.output_dir,
                                     replay=args.replay), indent=2, allow_nan=False))
    else:
        parser.error("provide trigger_json and output_dir, or --probe DIR")
