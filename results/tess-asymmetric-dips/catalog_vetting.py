"""Query the four catalog services for TIC positions from the frozen product table.

python catalog_vetting.py products.csv triggers.csv output_directory
Both CSVs must have a tic column; products.csv must also contain ra and dec (deg).
Use the corresponding product and trigger tables for each cohort.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


def query_url(base, params):
    return base + "?" + urllib.parse.urlencode(params)


def queries(ra, dec):
    def cone(radius, cols):
        return (f"1=CONTAINS(POINT('ICRS',{cols}),"
                f"CIRCLE('ICRS',{ra},{dec},{radius / 3600:.12f}))")

    def tap(base, sql):
        return query_url(base, dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=sql))

    gaia_sql = ("SELECT g.source_id,g.ra,g.dec,g.ref_epoch,g.pmra,g.pmdec,g.parallax,"
                "g.parallax_error,g.phot_g_mean_mag,g.bp_rp,g.ruwe,g.phot_variable_flag,"
                "g.non_single_star,g.has_epoch_photometry,v.best_class_name,v.best_class_score "
                "FROM gaiadr3.gaia_source AS g LEFT OUTER JOIN gaiadr3.vari_classifier_result "
                "AS v ON g.source_id=v.source_id WHERE " + cone(63, "g.ra,g.dec"))
    simbad_sql = ("SELECT b.main_id,b.otype,b.ra,b.dec,b.sp_type FROM basic AS b WHERE "
                  + cone(30, "b.ra,b.dec"))
    vsx = query_url("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv", {
        "-source": "B/vsx/vsx", "-c": f"{ra} {dec}", "-c.rs": 63,
        "-out": "OID,Name,Type,Period,RAJ2000,DEJ2000", "-out.add": "_r", "-out.max": 1000})
    columns = ("designation,ra,dec,w1mpro,w1sigmpro,w1snr,w2mpro,w2sigmpro,w2snr,"
               "w3mpro,w3sigmpro,w3snr,w4mpro,w4sigmpro,w4snr,cc_flags,ph_qual,"
               "ext_flg,j_m_2mass,h_m_2mass,k_m_2mass")
    wise = query_url("https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query", {
        "catalog": "allwise_p3as_psd", "spatial": "Cone", "objstr": f"{ra} {dec}",
        "radius": 12, "radunits": "arcsec", "selcols": columns, "outfmt": 3})
    return dict(gaia=tap("https://gea.esac.esa.int/tap-server/tap/sync", gaia_sql),
                simbad=tap("https://simbad.harvard.edu/simbad/sim-tap/sync", simbad_sql),
                vsx=vsx, wise=wise)


def decode(label, text):
    if label in ("gaia", "simbad"):
        key = "source_id" if label == "gaia" else "main_id"
        if not text.lower().startswith(key + ","):
            raise ValueError("unexpected TAP CSV header")
        return list(csv.DictReader(io.StringIO(text)))
    if label == "wise":
        table = ET.fromstring(text)
        keys = [f.get("name") for f in table.findall(".//{*}FIELD")]
        if table.find(".//{*}TABLEDATA") is None or "designation" not in keys:
            raise ValueError("missing AllWISE table")
        return [dict(zip(keys, [td.text or "" for td in tr.findall("{*}TD")]))
                for tr in table.findall(".//{*}TABLEDATA/{*}TR")]
    if "#Name: B/vsx/vsx" not in text or "STDB_ERR" in text or "#ERROR" in text.upper():
        raise ValueError("invalid VizieR table or database connection error")
    rows = []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) == 7 and cols[1].strip().isdigit():
            rows.append(dict(zip(["sep", "oid", "name", "type", "period", "ra", "dec"],
                                 [c.strip() for c in cols])))
    return rows


def fetch(item):
    tic, label, url = item
    result = dict(tic=tic, catalog=label, url=url)
    for attempt in range(3):
        result["utc"] = datetime.now(timezone.utc).isoformat()
        try:
            with urllib.request.urlopen(url, timeout=45) as h:
                data = h.read()
                result["http"] = h.status
            text = data.decode("utf-8")
            result.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), response=text)
            rows = decode(label, text)
            result.update(ok=True, count=len(rows))
            result.pop("error", None)
            return result, rows
        except (ValueError, OSError, ET.ParseError) as exc:
            result.update(ok=False, error=str(exc))
            if attempt < 2:
                time.sleep(2 ** attempt)
    return result, []


def separation(row, ra, dec):
    a, d = math.radians(float(row["ra"])), math.radians(float(row["dec"]))
    a0, d0 = math.radians(ra), math.radians(dec)
    v = math.sin((d-d0)/2)**2 + math.cos(d)*math.cos(d0)*math.sin((a-a0)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(v, 1))))*3600


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("products", type=Path)
    p.add_argument("triggers", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    with args.products.open() as f:
        products = {int(r["tic"]): (float(r["ra"]), float(r["dec"])) for r in csv.DictReader(f)}
    with args.triggers.open() as f:
        tics = sorted({int(r["tic"]) for r in csv.DictReader(f)})
    tasks = [(t, label, url) for t in tics for label, url in queries(*products[t]).items()]
    args.output.mkdir(parents=True, exist_ok=True)
    results = {t: {"tic": t} for t in tics}
    failures = Counter()
    with gzip.open(args.output / "catalog-receipts.jsonl.gz", "wt", encoding="utf-8") as out:
        with ThreadPoolExecutor(4) as pool:
            for receipt, rows in pool.map(fetch, tasks):
                out.write(json.dumps(receipt, separators=(",", ":")) + "\n")
                tic, label = receipt["tic"], receipt["catalog"]
                dest = results[tic]
                dest[label + "_status"] = "ok" if receipt["ok"] else "unavailable"
                if not receipt["ok"]:
                    failures[label] += 1
                    continue
                dest[label + "_count"] = len(rows)
                if label == "vsx":
                    dest["vsx_matches"] = json.dumps(rows, separators=(",", ":"))
                elif rows:
                    rows.sort(key=lambda r: separation(r, *products[tic]))
                    closest = rows[0]
                    dest[label + "_closest_arcsec"] = round(separation(closest, *products[tic]), 3)
                    if label == "gaia":
                        for k in ("source_id", "phot_g_mean_mag", "bp_rp", "parallax", "ruwe",
                                  "phot_variable_flag", "non_single_star", "has_epoch_photometry",
                                  "best_class_name"):
                            dest[k] = closest[k]
                        peers = [r for r in rows[1:] if r["phot_g_mean_mag"]]
                        if peers and closest["phot_g_mean_mag"]:
                            bright = min(peers, key=lambda r: float(r["phot_g_mean_mag"]))
                            dest["brightest_neighbor_source"] = bright["source_id"]
                            dest["brightest_neighbor_delta_G"] = round(
                                float(bright["phot_g_mean_mag"]) - float(closest["phot_g_mean_mag"]), 3)
                            dest["brightest_neighbor_sep_arcsec"] = round(separation(bright, *products[tic]), 3)
                    if label == "simbad":
                        for k in ("main_id", "otype", "sp_type"):
                            dest[k] = closest[k]
                    if label == "wise":
                        for k in ("designation", "ph_qual", "cc_flags", "ext_flg", "w1mpro",
                                  "w2mpro", "w3mpro", "w4mpro", "w3snr", "w4snr", "k_m_2mass"):
                            dest[k] = closest[k]
    with (args.output / "catalog-summary.csv").open("w", newline="") as out:
        fields = list(dict.fromkeys(k for row in results.values() for k in row))
        writer = csv.DictWriter(out, fields)
        writer.writeheader()
        writer.writerows(results.values())
    if failures:
        raise SystemExit(f"Incomplete catalog responses: {dict(failures)}")


if __name__ == "__main__":
    main()
