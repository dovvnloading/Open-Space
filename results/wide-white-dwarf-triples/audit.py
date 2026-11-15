"""Audit published Gaia identifiers against the cited release. Requires Python 3.9+."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
PUBLISHER = "https://content.cld.iop.org/journals/0004-637X/934/2/148/revision1/apjac78d9t2_mrt.txt"
SOURCE_SHA = "a741e22aac89e2bdf28c5376a44e380cc48b7de2bac686b5e983a93f74dc9128"
TAP = "https://gaia.ari.uni-heidelberg.de/tap/sync"
FIELDS = ["pair_row", "system", "component", "published_id", "ra", "dec", "g",
          "matched_id", "gaia_ra", "gaia_dec", "gaia_g", "cone_count"]
MATCH_RADIUS = 0.54  # arcsec, at the common reference epoch J2016.0


def rows_from_publisher(data):
    if hashlib.sha256(data).hexdigest() != SOURCE_SHA:
        raise ValueError("Publisher table checksum differs from the audited release")
    for number, line in enumerate((s for s in data.decode().splitlines()
                                    if s.startswith("WDWD")), 1):
        for component, offset in [("P", 0), ("S", 604)]:
            yield dict(pair_row=str(number), system=line[:36].strip(), component=component,
                       published_id=line[37+offset:56+offset].strip(),
                       ra=line[57+offset:68+offset].strip(),
                       dec=line[69+offset:80+offset].strip(),
                       g=line[125+offset:132+offset].strip())


def separation(row, source):
    # Local great-circle haversine formula, degrees to arcseconds.
    ra, dec, a, d = (math.radians(float(v)) for v in
                     (row["ra"], row["dec"], source["ra"], source["dec"]))
    h = math.sin((d-dec)/2)**2 + math.cos(dec)*math.cos(d)*math.sin((a-ra)/2)**2
    return math.degrees(2*math.asin(min(1, math.sqrt(h))))*3600


def fetch_gaia(condition):
    query = ("SELECT source_id,ra,dec,phot_g_mean_mag FROM gaiaedr3.gaia_source WHERE "
             + condition)
    url = TAP + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    with urlopen(url, timeout=120) as reply:
        reader = csv.DictReader(io.StringIO(reply.read().decode()))
        if reader.fieldnames != ["source_id", "ra", "dec", "phot_g_mean_mag"]:
            raise ValueError("Unexpected TAP response")
        return list(reader)


def cone_lookup(rows):
    sources = {}
    for i in range(0, len(rows), 8):
        cones = [f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r['ra']},{r['dec']},{MATCH_RADIUS/3600}))"
                 for r in rows[i:i+8]]
        for source in fetch_gaia(" OR ".join(cones)):
            sources[source["source_id"]] = source
    for row in rows:
        matches = [s for s in sources.values() if separation(row, s) <= MATCH_RADIUS]
        if len(matches) != 1:
            raise ValueError(f"Ambiguous or missing counterpart {row}: {len(matches)}")
        source = matches[0]
        yield dict(**row, matched_id=source["source_id"], gaia_ra=source["ra"],
                   gaia_dec=source["dec"], gaia_g=source["phot_g_mean_mag"], cone_count="1")


def binary_controls(rows):
    ids = sorted({r["published_id"] for r in rows})
    sources = {}
    for i in range(0, len(ids), 120):
        condition = "source_id IN (" + ",".join(ids[i:i+120]) + ")"
        for s in fetch_gaia(condition):
            sources[s["source_id"]] = s
    return control_summary(rows, sources)


def control_summary(rows, sources):
    offsets = [separation(r, sources[r["published_id"]]) for r in rows]
    mags = [abs(float(r["g"])-float(sources[r["published_id"]]["phot_g_mean_mag"])) for r in rows]
    return dict(components=len(rows), max_offset_arcsec=round(max(offsets), 9),
                max_abs_g_difference_mag=round(max(mags), 7),
                inconsistent_positions=sum(d > MATCH_RADIUS for d in offsets))


def summary(rows):
    mismatches = [r for r in rows if r["matched_id"] != r["published_id"]]
    max_offset = max(separation(r, {"ra":r["gaia_ra"], "dec":r["gaia_dec"]}) for r in rows)
    max_mag = max(abs(float(r["g"])-float(r["gaia_g"])) for r in rows)
    return dict(pair_rows=len({r["pair_row"] for r in rows}), component_records=len(rows),
                unique_stars=len({r["matched_id"] for r in rows}),
                mismatched_fields=len(mismatches), affected_rows=len({r["pair_row"] for r in mismatches}),
                affected_systems=len({r["system"] for r in mismatches}),
                max_offset_arcsec=round(max_offset, 9), max_abs_g_difference_mag=round(max_mag, 7))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--online", action="store_true", help="Requery the full publisher table and Gaia EDR3")
    args = parser.parse_args()
    with (HERE/"triple-components.csv").open() as fp:
        frozen = list(csv.DictReader(fp))
    if args.online:
        with urlopen(PUBLISHER, timeout=120) as reply:
            rows = list(rows_from_publisher(reply.read()))
        triples = [r for r in rows if r["system"].split()[0] != "WDWD"]
        new = list(cone_lookup(triples))
        if new != frozen:
            raise ValueError("Gaia or publisher rows differ from frozen triple projections")
        print(json.dumps({"binary_control":binary_controls([r for r in rows if r["system"].startswith("WDWD ")])}))
    print(json.dumps(summary(frozen), indent=2))


if __name__ == "__main__":
    main()
