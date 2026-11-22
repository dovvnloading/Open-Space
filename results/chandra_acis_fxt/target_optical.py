"""Targeted, post-search optical/IR counterpart evidence for two frozen X-ray flares."""
import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import xml.etree.ElementTree as ET

from astropy.time import Time
from archive import dump
from vet import CAMPAIGN_REF, GAIA, MIRRORS, decode_table, number, pm_position, read_input, separation

# These identify the real/0 rows in Actions run 44985879443, worker shards 31 and 19.
TARGETS = {
    23324: ("507a5aaf26b492ef2b5120242d9d33e281803667c8c1edc125fa905e6a09bbd0",
            127.76903026837093, 9.70538430082033, "599875677054390784"),
    32380: ("3c5f46567187114d7deb741d2a709bbc4e3566801a90a9cdb3f23a7c66a944cf",
            243.548594555635, -14.359688117051375, None),
}
CONE = 15.0  # All neighbors retained, in arcseconds, for association and blending review.
ROW_CAP = 1000
PS1 = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2"
IRSA = "https://irsa.ipac.caltech.edu/TAP/sync"
GAIA_FIELDS = """source_id ra dec ra_error dec_error ref_epoch astrometric_params_solved parallax parallax_error
pmra pmra_error pmdec pmdec_error ruwe duplicated_source astrometric_excess_noise
phot_g_mean_mag phot_bp_mean_mag phot_rp_mean_mag phot_g_mean_flux_over_error phot_bp_mean_flux_over_error
phot_rp_mean_flux_over_error bp_rp phot_bp_rp_excess_factor phot_variable_flag""".split()
WISE_FIELDS = """designation ra dec sigra sigdec w1mpro w1sigmpro w2mpro w2sigmpro w3mpro w3sigmpro
w4mpro w4sigmpro w1snr w2snr w3snr w4snr ph_qual cc_flags ext_flg nb na var_flg""".split()
DOCS = {"gaia": "https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html",
        "panstarrs": "https://catalogs.mast.stsci.edu/docs/panstarrs.html",
        "allwise": "https://irsa.ipac.caltech.edu/data/WISE/docs/release/AllWISE/expsup/sec2_1a.html"}


def diff(values, a, b, error_a=None, error_b=None):
    x, y = number(values.get(a)), number(values.get(b))
    if x is None or y is None or min(x, y) <= -900:
        return None
    ea, eb = number(values.get(error_a)), number(values.get(error_b))
    return {"value_mag": x-y, "sigma_independent_mag": math.hypot(ea, eb)
            if ea is not None and eb is not None and min(ea, eb) > 0 else None}


def tap_result(label, endpoint, query, fields, src, replay, receipts):
    raw = read_input(label, endpoint, src, replay, receipts, {
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable", "MAXREC": ROW_CAP, "QUERY": query})
    rows = decode_table(raw, fields, None)
    if len(rows) >= ROW_CAP:
        raise ValueError("TAP result reached the requested row ceiling")
    meta = {}
    for f in ET.fromstring(raw).iter():
        if f.tag.rsplit("}", 1)[-1] == "FIELD":
            meta[f.get("name", "").lower()] = dict(f.attrib)
    for c in ("ra", "dec"):
        if str(meta.get(c, {}).get("unit", "")).lower() not in {"deg", "degree", "degrees"}:
            raise ValueError(f"Missing degree unit on TAP column {c}")
    return rows, meta


def tap_service(name, lead, src, replay, out, known=None, prior_gaia_id=None):
    ra, dec = lead["ra_deg"], lead["dec_deg"]
    gaia = name.startswith("gaia")
    fields = GAIA_FIELDS if gaia else WISE_FIELDS
    table = "gaiadr3.gaia_source" if gaia else "allwise_p3as_psd"
    predicate = (f"source_id={known}" if known else
                 f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.12f},{dec:.12f},{CONE/3600:.12f}))")
    query = f"SELECT {','.join(fields)} FROM {table} WHERE {predicate}"
    endpoints = (GAIA, MIRRORS["gaia"][1]) if gaia else (IRSA,)
    out.update(table=table, query=query, errors=[], catalog_docs=DOCS["gaia" if gaia else "allwise"])
    for i, endpoint in enumerate(endpoints):
        try:
            rows, meta = tap_result(f"{name}_{i}", endpoint, query, fields, src, replay, out["receipts"])
        except Exception as error:
            out["errors"].append({"endpoint": endpoint, "error": str(error)})
            continue
        if known and (len(rows) != 1 or str(rows[0]["source_id"]) != known):
            out["errors"].append({"endpoint": endpoint, "error": "Expected Gaia ID absent or nonunique"})
            continue
        out.update(columns=meta, endpoint=endpoint)
        for row in rows:
            if gaia:
                x, y = pm_position(row, Time(lead["event_mjd"], format="mjd", scale="utc").jyear)
                derived = {"matches_prior_gaia_id": str(row["source_id"]) == prior_gaia_id if prior_gaia_id else None,
                           "offset_catalog_epoch_arcsec": separation(ra, dec, float(row["ra"]), float(row["dec"])),
                           "motion_propagated": all(number(row[x]) is not None for x in ("pmra", "pmdec")),
                           "offset_event_epoch_arcsec": separation(ra, dec, x, y)}
            else:
                x, y = float(row["ra"]), float(row["dec"])
                quality = row["ph_qual"] or ""
                color = (diff(row, "w1mpro", "w2mpro", "w1sigmpro", "w2sigmpro")
                         if len(quality) == 4 and all(f in "ABC" for f in quality[:2]) else None)
                derived = {"w1_minus_w2_vega": color}
            offset = separation(ra, dec, x, y)
            out["rows"].append({"offset_arcsec": offset, "within_xray_r90": offset <= lead["position_error_arcsec"],
                                "derived": derived, "catalog_values": row})
        return
    raise ValueError("No endpoint returned the required TAP result")


def ps1_service(kind, lead, src, replay, out):
    url = f"{PS1}/{kind}"
    out.update(catalog_docs=DOCS["panstarrs"], table=f"PS1 DR2 {kind}")
    metadata = json.loads(read_input(f"ps1_{kind}_schema", url+"/metadata.json", src, replay, out["receipts"]))
    columns = {c["name"].lower(): c for c in metadata}
    suffix = "mean" if kind == "mean" else "stack"
    rac, dec = "ra"+suffix, "dec"+suffix
    prefix = "mean" if kind == "mean" else ""
    bands = "grizy"
    mags = [b+prefix+c for b in bands for c in ("psfmag", "psfmagerr", "kronmag", "kronmagerr")]
    required = ["objid", rac, dec, "qualityflag", "objinfoflag", *mags]
    if kind == "stack":
        required += ["primarydetection", "bestdetection"]
    missing = set(required)-columns.keys()
    if missing:
        raise ValueError(f"Required PS1 fields absent: {sorted(missing)}")
    optional = [rac+"err", dec+"err", "epochmean", "ndetections",
                *("n"+b for b in bands), *(b+q for b in bands for q in
                (("flags", "qfperfect", "meanpsfmagnpt", "meankronmagnpt") if kind == "mean" else
                 ("infoflag", "infoflag2", "infoflag3", "psfqf", "psfmajorfwhm", "psfminorfwhm")))]
    selected = required + [c for c in optional if c in columns]
    out.update(columns={c: columns[c] for c in selected}, missing_optional=sorted(set(optional)-columns.keys()))
    params = {"ra": lead["ra_deg"], "dec": lead["dec_deg"], "radius": CONE/3600, "pagesize": ROW_CAP,
              "columns": "["+",".join(columns[c]["name"] for c in selected)+"]"}
    raw = read_input(f"ps1_{kind}", url+".csv", src, replay, out["receipts"], params)
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    headers = [k.lower() for k in reader.fieldnames or []]
    if set(headers) != set(selected) or len(headers) != len(selected):
        raise ValueError("PS1 CSV has missing, duplicate or unexpected columns")
    for row in reader:
        if None in row or None in row.values():
            raise ValueError("Malformed PS1 CSV row")
        vals = {k.lower(): v for k, v in row.items()}
        coords = [number(vals[c]) for c in (rac, dec)]
        offset = (separation(lead["ra_deg"], lead["dec_deg"], *coords)
                  if all(v is not None and v > -900 for v in coords) else None)
        color = {}
        for a, b in zip(bands, bands[1:]):
            m, n = a+prefix+"psfmag", b+prefix+"psfmag"
            color[a+"-"+b] = diff(vals, m, n, m+"err", n+"err")
        shape = {}
        for band in bands:
            m, n = band+prefix+"psfmag", band+prefix+"kronmag"
            shape[band] = diff(vals, m, n, m+"err", n+"err")
        out["rows"].append({"catalog_values": vals, "offset_arcsec": offset,
                            "within_xray_r90": None if offset is None else offset <= lead["position_error_arcsec"],
                            "psf_color_ab": color, "psf_minus_kron_mag": shape})
    if len(out["rows"]) >= ROW_CAP:
        raise ValueError("PS1 result reached page boundary")


def main(args):
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != CAMPAIGN_REF or
            os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("Use the Chandra research branch Actions workflow")
    body = args.leads.read_bytes()
    src = args.out / "source"
    src.mkdir(parents=True, exist_ok=True)
    summary = {"leads_sha256": hashlib.sha256(body).hexdigest(), "cone_arcsec": CONE, "leads": [],
               "status": "processing", "expected_obsids": sorted(TARGETS),
               "limitations": "Ancillary targeted follow-up. Positive source rows and offsets are descriptive. "
               "The X-ray r90 is approximate; PS1 and WISE coordinates remain at catalog epoch. Empty cones "
               "do not give image-based limiting magnitudes or prove image coverage. PS1 mean/stack can repeat "
               "objects. Missing and -999 magnitudes are excluded from differences. PS1 colors use AB, WISE "
               "uses Vega; no extinction correction. Color errors assume independent measurements. PSF and "
               "Kron errors are correlated, so their quadrature is only a reference scale. Quality, blending "
               "and saturation flags remain in the raw rows. These data alone do not identify source nature."}
    leads = json.loads(body)
    if not isinstance(leads, list):
        raise ValueError("Expected a list of normalized lead records")
    for lead in leads:
        oid = int(lead["obsid"])
        if oid not in TARGETS or lead["trigger_key"] != "real/0":
            continue
        sha, ra, dec, known = TARGETS[oid]
        if lead["input_sha256"] != sha or (lead["ra_deg"], lead["dec_deg"]) != (ra, dec):
            raise ValueError(f"ObsID {oid} does not match frozen event product and centroid")
        if any(r["obsid"] == oid for r in summary["leads"]):
            raise ValueError(f"Duplicate normalized lead for {oid}")
        record = {k: lead[k] for k in ("obsid", "trigger_key", "input_sha256", "ra_deg", "dec_deg",
                                       "position_error_arcsec", "event_mjd")}
        record["obsid"] = oid
        record["catalogs"] = {}
        summary["leads"].append(record)
        for name in ["gaia_cone", *(["gaia_prior_id"] if known else []), "ps1_mean", "ps1_stack", "allwise"]:
            out = {"status": "incomplete", "rows": [], "receipts": {}}
            record["catalogs"][name] = out
            try:
                if name.startswith("ps1_"):
                    ps1_service(name[4:], lead, src, args.replay, out)
                else:
                    tap_service(name, lead, src, args.replay, out, known if name == "gaia_prior_id" else None, known)
                out["rows"].sort(key=lambda x: x["offset_arcsec"] if x["offset_arcsec"] is not None else math.inf)
                out["status"] = "rows_returned" if out["rows"] else "empty_validated_cone"
            except Exception as error:
                out["error"] = str(error)
            dump(args.out / "review.json", summary)
            print(f"{oid} {name} {out['status']} ({len(out['rows'])} rows)", flush=True)
    summary["missing_obsids"] = sorted(set(TARGETS)-{r["obsid"] for r in summary["leads"]})
    summary["status"] = ("incomplete" if summary["missing_obsids"] or any(
        c["status"] == "incomplete" for r in summary["leads"] for c in r["catalogs"].values()) else "complete")
    dump(args.out / "review.json", summary)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--leads", type=Path, required=True)
    cli.add_argument("--out", type=Path, required=True)
    cli.add_argument("--replay", action="store_true")
    main(cli.parse_args())
