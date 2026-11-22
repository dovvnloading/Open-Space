"""Describe Gaia DR3 sources already associated with real Chandra leads by lead_optical."""
import argparse
from collections import defaultdict
import json
import math
import os
from pathlib import Path
import re
from urllib.parse import urlencode
from xml.etree import ElementTree as ET

from astropy.time import Time

import campaign
import ingest_primary
from target_optical import DOCS
from vet import GAIA, MIRRORS, decode_table, digest, pm_position, read_input, save_bytes, separation


ENDPOINTS = (GAIA, MIRRORS["gaia"][1])
SOURCE = "gaiadr3.gaia_source"
VAR = "gaiadr3.vari_classifier_result"
BATCH = 40
CAP = 2000
DSC = tuple("classprob_dsc_combmod_" + c for c in ("star", "galaxy", "quasar", "whitedwarf", "binarystar"))
PARAMS = tuple(p + "_gspphot" + s for p in ("teff", "logg", "mh", "distance", "azero", "ag", "ebpminrp")
               for s in ("", "_lower", "_upper"))
COLS = {SOURCE: ("source_id", "solution_id", "phot_g_mean_mag", "bp_rp", "parallax", "parallax_error",
                 "pmra", "pmra_error", "pmdec", "pmdec_error", "astrometric_params_solved", "ruwe",
                 "visibility_periods_used", "astrometric_excess_noise", "duplicated_source", "phot_variable_flag", "has_epoch_photometry", "has_xp_continuous",
                 "has_xp_sampled", "has_rvs", "has_mcmc_gspphot", "in_qso_candidates", "in_galaxy_candidates",
                 "flags_gspphot", "libname_gspphot", *DSC, *PARAMS),
        VAR: ("source_id", "solution_id", "classifier_name", "best_class_name", "best_class_score")}
SCHEMA_FIELDS = ("table_name", "column_name", "datatype", "unit", "description")
SCHEMA_QUERY = ("SELECT " + ",".join(SCHEMA_FIELDS) + " FROM TAP_SCHEMA.columns WHERE table_name IN "
                "('gaiadr3.gaia_source','gaiadr3.vari_classifier_result')")
VAR_DOC = ("https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/"
           "sec_dm_variability_tables/ssec_dm_vari_classifier_result.html")


def sid(value):
    if type(value) not in (int, str) or not re.fullmatch(r"[1-9][0-9]{0,18}", str(value)) or int(value) >= 2**63:
        raise ValueError("Expected an exact positive signed-64-bit Gaia source ID")
    return str(value)


def associations(root):
    base = [root/"lead_optical"/"observations"]
    base.extend((root/"pointed_visits").glob("**/optical/observations"))
    leads = {}
    for directory in sorted(set(base)):
        for path in sorted(directory.glob("*/*/real/*.json.gz")):
            obj = campaign.read_obs(path)
            lead = obj["lead"]
            oid, h, tag = (lead[k] for k in ("obsid", "input_sha256", "trigger_key"))
            if (obj["schema"] != 1 or type(oid) is not int or not re.fullmatch(r"[a-f0-9]{64}", h)
                    or not re.fullmatch(r"real/(0|[1-9][0-9]*)", tag)
                    or path != directory/str(oid)/h/(tag+".json.gz")):
                raise ValueError(f"Invalid optical record identity {path}")
            rec = leads.setdefault((oid, h, tag), {"lead": lead, "optical_records": [], "associations": {}})
            if rec["lead"] != lead:
                raise ValueError(f"Conflicting Chandra lead coordinates {path}")
            rel = path.relative_to(root).as_posix()
            prior = {sid(x) for x in obj["prior_gaia_ids"]}
            cats = {k: v for k, v in obj["catalogs"].items() if k == "gaia_cone" or k.startswith("gaia_id_")}
            if set(cats) != {"gaia_cone", *("gaia_id_"+s for s in prior)}:
                raise ValueError(f"Unexpected optical service list {path}")
            rec["optical_records"].append({"path": rel, "sha256": digest(path.read_bytes()),
                "complete": obj["complete"], "gaia_service_status": {k: v["status"] for k, v in cats.items()}})
            ids = prior | {sid(row["catalog_values"]["source_id"]) for cat in cats.values() for row in cat["rows"]
                           if row["within_xray_r90"] is True}
            for s in sorted(ids):
                assoc = rec["associations"].setdefault(s, {"primary_match_records": [], "positions": []})
                if s in prior:
                    assoc["primary_match_records"].append(rel)
                for name, cat in cats.items():
                    for row in cat["rows"]:
                        values = row["catalog_values"]
                        if sid(values["source_id"]) != s:
                            continue
                        xy = pm_position(values, Time(lead["event_mjd"], format="mjd", scale="utc").jyear)
                        off = separation(lead["ra_deg"], lead["dec_deg"], *xy)
                        if (not math.isclose(off, row["offset_arcsec"], abs_tol=1e-7) or
                                (off <= lead["position_error_arcsec"]) != row["within_xray_r90"]):
                            raise ValueError(f"Inconsistent stored Gaia position or aperture {path}")
                        assoc["positions"].append({"record": rel, "catalog": name, "offset_arcsec": row["offset_arcsec"],
                            "within_xray_r90": row["within_xray_r90"], "epoch_ra_deg": xy[0], "epoch_dec_deg": xy[1],
                            "ref_epoch": values["ref_epoch"], "catalog_ra_deg": values["ra"], "catalog_dec_deg": values["dec"],
                            "g_mag": values.get("phot_g_mean_mag"), **row["derived"]})
    return [leads[k] for k in sorted(leads)]


def parameters(query):
    return {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "MAXREC": CAP, "QUERY": query}


def query_ids(table, fields, ids):
    if table not in COLS or fields != [c for c in COLS[table] if c in fields] or not 0 < len(ids) <= BATCH:
        raise ValueError("Invalid batch contract")
    return f"SELECT {','.join(fields)} FROM {table} WHERE source_id IN ({','.join(sid(x) for x in ids)})"


def schema_columns(rows, table):
    subset = [r for r in rows if r["table_name"].lower() == table]
    meta = {r["column_name"].lower(): r for r in subset}
    required = {"source_id", "phot_g_mean_mag"} if table == SOURCE else {"source_id", "best_class_name", "best_class_score"}
    if len(subset) != len(meta) or not required <= meta.keys():
        raise ValueError(f"{table} required columns unavailable or duplicated")
    datatype = meta["source_id"]["datatype"].upper().removeprefix("ADQL:")
    if datatype not in {"LONG", "BIGINT", "INT8", "INTEGER64"}:
        raise ValueError("Unverified 64-bit source_id schema")
    return [c for c in COLS[table] if c in meta]


def read_table(raw, fields, ids=None, table=None):
    fields_xml = [e for e in ET.fromstring(raw).iter() if e.tag.rsplit("}", 1)[-1] == "FIELD"]
    names = [e.get("name", "").lower() for e in fields_xml]
    if len(names) != len(set(names)) or set(names) != set(fields):
        raise ValueError("Unexpected TAP result columns")
    if ids is not None and fields_xml[names.index("source_id")].get("datatype", "").lower() != "long":
        raise ValueError("Gaia source_id must be encoded as a VOTable long")
    rows = decode_table(raw, fields, None)
    if len(rows) >= CAP:
        raise ValueError("Gaia result reached its row ceiling")
    if ids is None:
        return rows
    grouped = {s: [] for s in ids}
    for row in rows:
        s = sid(row["source_id"])
        if s not in grouped or (table == SOURCE and grouped[s]):
            raise ValueError("Unrequested or repeated Gaia source row")
        for k, v in row.items():
            if v is None or v == "" or v.lower() == "nan":
                row[k] = None
            elif k in DSC or k == "best_class_score":
                if not 0 <= float(v) <= 1:
                    raise ValueError("Invalid Gaia probability or classification score")
        if table == VAR and any(r.get("classifier_name") == row.get("classifier_name") for r in grouped[s]):
            raise ValueError("Repeated Gaia variability-classifier row")
        grouped[s].append(row)
    return grouped


def body(out, receipt):
    """Validate the request contract and original HTTP bytes, including failed attempts."""
    if receipt["endpoint"] not in ENDPOINTS or not re.fullmatch(r"[0-9]+-[0-9]+", receipt["run_attempt"]):
        raise ValueError("Unexpected recorded endpoint or execution identity")
    memo = receipt["http"]
    if memo["url"] != receipt["endpoint"] + "?" + urlencode(parameters(receipt["query"])):
        raise ValueError("Changed Gaia query receipt")
    raw = None
    for a in memo["attempts"]:
        name = a.get("body_file")
        if not name:
            continue
        if Path(name).name != name:
            raise ValueError("Unsafe Gaia response filename")
        data = (out/"source"/receipt["run_attempt"]/name).read_bytes()
        if len(data) != a["size_bytes"] or digest(data) != a["sha256"]:
            raise ValueError("Altered Gaia response")
        if name == memo.get("body_file"):
            raw = data
    if memo.get("ok") and (raw is None or digest(raw) != memo["sha256"]):
        raise ValueError("Missing successful Gaia response")
    return raw


def fetch(out, endpoint, query, fields, replay, run, *, ids=None, table=None):
    contract = {"query": query, "endpoint": endpoint, "table": table, "ids": ids, "fields": list(fields)}
    key = digest(json.dumps(contract, sort_keys=True).encode())
    home = out/"requests"/key
    ok = home/"ok.json"
    if ok.is_file():
        rec = campaign.read(ok)
        if any(rec[k] != v for k, v in contract.items()):
            raise ValueError("Changed Gaia cache contract")
        return read_table(body(out, rec), fields, ids, table), ok.relative_to(out).as_posix(), None
    if replay:
        result = None
        for p in sorted(home.glob("error-*.json")):
            failed = campaign.read(p)
            if any(failed[k] != v for k, v in contract.items()):
                raise ValueError("Changed failed Gaia contract")
            if failed.get("http"):
                body(out, failed)
            result = p.relative_to(out).as_posix()
        return None, result, "No validated response"
    home.mkdir(parents=True, exist_ok=True)
    src = out/"source"/run
    src.mkdir(parents=True, exist_ok=True)
    http = {}
    rec = {**contract, "run_attempt": run}
    try:
        raw = read_input("gaia", endpoint, src, False, http, parameters(query))
        rows = read_table(raw, fields, ids, table)
    except (ValueError, OSError, ET.ParseError, TypeError, IndexError, KeyError) as exc:
        rec.update(http=http.get("gaia"), error=f"{type(exc).__name__}: {exc}")
        p = home/f"error-{run}.json"
        campaign.dump(p, rec)
        return None, p.relative_to(out).as_posix(), "No validated response"
    rec["http"] = http["gaia"]
    campaign.dump(ok, rec)
    return rows, ok.relative_to(out).as_posix(), None


def lookup(out, ids, replay, run):
    sources = {s: {} for s in ids}
    schemas = []
    # Prefer the ESA responses. Use AIP's DR3 mirror only for still-unanswered IDs.
    for endpoint in ENDPOINTS if ids else ():
        if all(len(v) == len(COLS) for v in sources.values()):
            break
        schema, ref, error = fetch(out, endpoint, SCHEMA_QUERY, SCHEMA_FIELDS, replay, run)
        description = {"endpoint": endpoint, "schema_receipt": ref, "error": error, "tables": {}}
        schemas.append(description)
        if schema is None:
            continue
        for table in COLS:
            needed = {s for s in ids if table not in sources[s]}
            if not needed:
                continue
            meta = {"failures": []}
            description["tables"][table] = meta
            try:
                fields = schema_columns(schema, table)
            except ValueError as exc:
                meta["schema_error"] = str(exc)
                continue
            meta["selected"] = fields
            meta["missing_optional"] = [c for c in COLS[table] if c not in fields]

            def apply(batch):
                rows, pointer, err = fetch(out, endpoint, query_ids(table, fields, batch), fields,
                                           replay, run, ids=batch, table=table)
                if rows is None:
                    meta["failures"].append({"ids": batch, "receipt": pointer, "error": err})
                    return
                for s in needed.intersection(rows):
                    status = "row" if rows[s] else ("no_returned_source_row" if table == SOURCE else "no_published_row")
                    sources[s][table] = {"status": status, "rows": rows[s],
                                         "receipt": pointer, "schema_receipt": ref}
                needed.difference_update(rows)

            for record in sorted((out/"requests").glob("*/ok.json")):
                previous = campaign.read(record)
                batch = previous["ids"]
                if (previous["endpoint"] == endpoint and previous["table"] == table and previous["fields"] == fields
                        and batch and needed.intersection(batch)):
                    apply(batch)
            remaining = sorted(needed, key=int)
            for start in range(0, len(remaining), BATCH):
                apply(remaining[start:start+BATCH])
    for result in sources.values():
        for table in COLS:
            result.setdefault(table, {"status": "unavailable", "rows": None, "receipt": None, "schema_receipt": None})
    return sources, schemas


def markdown(report, root, out):
    sources = report["sources"]
    leads = report["leads"]
    n_assoc = sum(len(x["associations"]) for x in leads)
    missing = sum(v["status"] == "unavailable" for s in sources.values() for v in s.values())
    lines = ["# Gaia DR3 properties of optical sources at Chandra leads", "",
             f"{len(leads)} deposited Chandra leads have {n_assoc} source associations with "
             f"{len(sources)} distinct Gaia DR3 identifiers. {missing} source-table lookups remain unavailable.", "",
             "Each association comes from the primary Gaia match or an optical position inside the X-ray "
             "90% error circle at the event epoch. Positions without both proper-motion components use the "
             "Gaia catalog coordinates. A positional coincidence alone cannot establish the origin of a flare.", "",
             "DSC values are Gaia posterior class probabilities and include the global class prior. "
             "GSP-Phot assumes a single star. Brackets give the published 16th and 84th percentiles. "
             "The variability score is a normalized rank. A missing variability row leaves variability unconstrained. "
             "These descriptions do not enter either the blind selection or the held-out test. "
             f"[Gaia source data model]({DOCS['gaia']}) · [Variability table]({VAR_DOC})", "",
             "Each lookup links to its exact request and response checksum. The association ledger records "
             "every event-epoch position, proper-motion flag, input hash, ancillary flag, physical parameter, "
             "and schema or transport failure. An em dash marks an unavailable or unpublished value. "
             "[Full association ledger](associations.json)", "",
             "| Chandra lead | Gaia DR3 source | Offset (arcsec) | G (mag) | P(star / galaxy / quasar) | "
             "T eff (K) | Distance (pc) | Variability class and score | Source / variability lookup |",
             "|---|---|---:|---:|---|---|---|---|---|"]
    def show(v):
        return "—" if v is None else str(v).replace("|", "&#124;")

    def interval(row, key):
        median, a, b = (row.get(key+s) for s in ("", "_lower", "_upper"))
        return "—" if all(v is None for v in (median, a, b)) else f"{show(median)} [{show(a)}, {show(b)}]"

    def lookup_link(rec):
        label = {"row": "reported", "no_published_row": "no published row",
                 "no_returned_source_row": "expected source absent", "unavailable": "unavailable"}[rec["status"]]
        if not rec["receipt"]:
            return label
        packet = campaign.read(out/rec["receipt"])
        raw = (Path("source")/packet["run_attempt"]/packet["http"]["body_file"]).as_posix()
        return f"[{label}]({rec['receipt']}) ([bytes]({raw}))"

    for item in leads:
        lead = item["lead"]
        link = os.path.relpath(root/item["optical_records"][0]["path"], out)
        name = f"[{lead['obsid']} {lead['trigger_key']}]({link})"
        if not item["associations"]:
            gaia_ok = any(all(s != "incomplete" for s in rec["gaia_service_status"].values())
                          for rec in item["optical_records"])
            lines.append(f"| {name} | {'No linked source' if gaia_ok else 'Optical lookup incomplete'} | — | — | — | — | — | — | — |")
        for s, assoc in sorted(item["associations"].items(), key=lambda pair: int(pair[0])):
            rec, var = sources[s][SOURCE], sources[s][VAR]
            row = (rec["rows"] or [{}])[0]
            offsets = sorted({f"{p['offset_arcsec']:.3f}" + (" inside" if p["within_xray_r90"] else " outside")
                              for p in assoc["positions"]})
            flags = ", ".join(show(v.get("best_class_name")) + " " + show(v.get("best_class_score"))
                              for v in (var["rows"] or []))
            g = row.get("phot_g_mean_mag")
            if g is None:
                g = next((p["g_mag"] for p in assoc["positions"] if p["g_mag"] is not None), None)
            cells = [name, s, ", ".join(offsets) or "—", show(g), " / ".join(show(row.get(c)) for c in DSC[:3]),
                     interval(row, "teff_gspphot"), interval(row, "distance_gspphot"), flags or "—",
                     lookup_link(rec) + " / " + lookup_link(var)]
            lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def run(root, out, replay=False):
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REPOSITORY") != ingest_primary.REPO
            or os.getenv("GITHUB_REF") != "refs/heads/" + ingest_primary.BRANCH
            or os.getenv("GITHUB_WORKFLOW_REF", "").split("@", 1)[0] != ingest_primary.REPO + "/" + ingest_primary.FLOW):
        raise RuntimeError("Run in the Chandra research workflow")
    run_id = f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    if not re.fullmatch(r"[0-9]+-[0-9]+", run_id):
        raise ValueError("Invalid Actions execution identity")
    leads = associations(root)
    out.mkdir(parents=True, exist_ok=True)
    sources, schemas = lookup(out, sorted({s for v in leads for s in v["associations"]}, key=int), replay, run_id)
    paths = {x["schema_receipt"] for x in schemas if x["schema_receipt"]}
    for provider in schemas:
        for meta in provider["tables"].values():
            paths.update(f["receipt"] for f in meta["failures"] if f["receipt"])
    paths.update(v["receipt"] for s in sources.values() for v in s.values() if v["receipt"])
    report = {"schema": 1, "implementation_sha256": {n: digest((Path(__file__).parent/n).read_bytes())
              for n in ("gaia_properties.py", "vet.py", "lead_optical.py", "target_optical.py")},
              "leads": leads, "sources": sources, "providers": schemas,
              "receipt_sha256": {p: digest((out/p).read_bytes()) for p in sorted(paths)}}
    summary = out/"associations.json"
    page = out/"README.md"
    rendered = markdown(report, root, out).encode()
    if replay:
        if report != campaign.read(summary) or rendered != page.read_bytes():
            raise ValueError("Gaia characterization differs on offline replay")
    else:
        campaign.dump(summary, report)
        save_bytes(page, rendered)
    print(f"Gaia DR3: {len(leads)} deposited leads, {len(sources)} linked source IDs", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--replay", action="store_true", help="Verify without HTTP")
    args = parser.parse_args()
    run(args.evidence, args.out or args.evidence/"gaia_dr3_properties", args.replay)
