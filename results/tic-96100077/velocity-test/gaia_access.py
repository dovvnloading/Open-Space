"""Resolve the target across Gaia releases and retrieve its advertised DR4 RV products."""

import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit
from urllib.request import urlopen
from xml.etree import ElementTree as ET

DR3_ID = "5946882743944130304"
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
INDEX = "https://gea.esac.esa.int/data-server/datalink/links"
PRODUCTS = {
    "EPOCH_PARAMETERS_RVS_SINGLE": "single.csv",
    "EPOCH_PARAMETERS_RVS_DOUBLE": "double.csv",
}


def get(url, receipt):
    with urlopen(url, timeout=90) as reply:
        blob = reply.read()
        receipt.append(
            {
                "url": url,
                "status": reply.status,
                "bytes": len(blob),
                "sha256": hashlib.sha256(blob).hexdigest(),
            }
        )
    return blob


def query(sql, receipt):
    url = (
        TAP
        + "?"
        + urlencode(
            {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": sql}
        )
    )
    data = get(url, receipt)
    if data.lstrip().startswith(b"<"):
        raise ValueError("Expected TAP CSV: " + data[:350].decode(errors="replace"))
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))


def datalinks(release, sid, receipt):
    blob = get(INDEX + "?" + urlencode({"ID": f"Gaia DR{release} {sid}"}), receipt)
    root = ET.fromstring(blob)
    for info in root.findall(".//{*}INFO"):
        if info.get("name") == "QUERY_STATUS" and info.get("value") in {
            "ERROR",
            "OVERFLOW",
        }:
            raise ValueError(info.text or "DataLink query failed")
    rows = []
    for table in root.findall(".//{*}TABLE"):
        names = [f.get("name") for f in table.findall("{*}FIELD")]
        for row in table.findall(".//{*}TR"):
            cells = dict(zip(names, (c.text or "" for c in row.findall("{*}TD"))))
            rows.append(
                {
                    k: cells[k]
                    for k in ("ID", "description", "access_url", "error_message")
                }
            )
    return rows


def run(folder, chosen):
    folder.mkdir(parents=True, exist_ok=True)
    requests = []
    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "dr3_id": DR3_ID,
        "requests": requests,
        "retrieved": {},
    }
    report["dr3_source"] = query(
        "SELECT source_id,ref_epoch,ra,dec,pmra,pmdec,phot_g_mean_mag,"
        "grvs_mag,has_rvs,has_epoch_rv FROM gaiadr3.gaia_source "
        f"WHERE source_id={DR3_ID}",
        requests,
    )[0]
    report["dr3_links"] = datalinks(3, DR3_ID, requests)
    found = query(
        "SELECT table_name FROM tap_schema.tables WHERE "
        "table_name='gaiadr4.dr3_neighbourhood'",
        requests,
    )
    report["dr4_match_table_rows"] = found

    if not found:
        report["status"] = "DR4 neighbourhood table unavailable on public TAP"
        return report

    matches = query(
        "SELECT source_id,dr3_source_id,angular_distance,magnitude_difference,"
        "proper_motion_propagation,probable_match FROM gaiadr4.dr3_neighbourhood "
        f"WHERE dr3_source_id={DR3_ID}",
        requests,
    )
    report["dr4_neighbourhood"] = matches
    # DR4 reference epoch J2017.5; the target moves ~12 mas since J2016.
    # The larger cone includes nearby sources if the catalog split this object.
    cone = query(
        "SELECT source_id,ra,dec,ref_epoch,phot_g_mean_mag FROM gaiadr4.gaia_source "
        "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
        "CIRCLE('ICRS',264.526817,-49.349218,0.001667))",
        requests,
    )
    report["dr4_six_arcsec_cone"] = cone
    ids = {row["source_id"] for row in cone}
    if chosen:
        if chosen not in ids:
            raise ValueError("Selected DR4 source ID is absent from the target cone")
        sid = chosen
    else:
        probable = [
            r["source_id"]
            for r in matches
            if r["probable_match"].lower() in ("true", "1", "t")
        ]
        if len(probable) != 1 or probable[0] not in ids or len(matches) != 1:
            report["status"] = (
                "Inspect neighbourhood and cone for alternate matches or a source split"
            )
            return report
        sid = probable[0]
    report["dr4_id"] = sid
    rows = datalinks(4, sid, requests)
    report["dr4_links"] = rows
    for row in rows:
        if not row["access_url"] or row["error_message"]:
            continue
        parts = urlsplit(row["access_url"])
        params = dict(parse_qsl(parts.query))
        product = params.get("RETRIEVAL_TYPE", "").upper()
        if product not in PRODUCTS:
            continue
        if product in report["retrieved"]:
            raise ValueError(
                f"Multiple advertised URLs for {product}; inspect the index"
            )
        if row["ID"] != f"Gaia DR4 {sid}":
            raise ValueError(
                f"DataLink identifier differs from the selected source: {row['ID']}"
            )
        params.update(FORMAT="CSV", DATA_STRUCTURE="DATAMODEL_GAIA", VALID_DATA="false")
        url = urlunsplit(parts._replace(query=urlencode(params)))
        blob = get(url, requests)
        reader = csv.DictReader(io.StringIO(blob.decode("utf-8-sig")))
        if "source_id" not in (reader.fieldnames or ()):
            raise ValueError(f"{product} response is not an epoch parameter CSV")
        n = 0
        for line in reader:
            if line["source_id"] != sid:
                raise ValueError("Mixed source IDs in downloaded product")
            n += 1
        path = folder / PRODUCTS[product]
        path.write_bytes(blob)
        report["retrieved"][product] = {
            "file": path.name,
            "rows": n,
            "sha256": requests[-1]["sha256"],
        }
    report["status"] = "DR4 index checked for selected source"
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, help="new or empty directory for this retrieval"
    )
    parser.add_argument(
        "--source-id", help="DR4 ID independently selected from the recorded cone"
    )
    args = parser.parse_args()
    if args.directory.exists() and any(args.directory.iterdir()):
        parser.error(
            "Use a new or empty directory to preserve prior receipts and source files"
        )
    result = run(args.directory, args.source_id)
    (args.directory / "access.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"])
    print(json.dumps(result["retrieved"]))
