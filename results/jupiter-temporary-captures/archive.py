"""Freeze or rebuild the JPL inputs. Keep raw response hashes with compact tables."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
VECTOR_WINDOWS = [
    ("uq146-past","2025 UQ146","1978-01-01","1993-01-01"),
    ("uq146-future","2025 UQ146","2062-01-01","2078-01-01"),
    ("oe31-past","2019 OE31","2005-01-01","2022-01-01"),
    ("oe31-future","2019 OE31","2087-01-01","2108-01-01"),
    ("positive-467p","467P","2002-01-01","2018-01-01"),
    ("negative-255p","255P","1899-01-01","1903-01-01"),
]


def requests():
    base = "https://ssd-api.jpl.nasa.gov/"
    cad = {"body":"Juptr","date-min":"1900-01-01","date-max":"2100-01-01",
           "dist-max":0.36,"neo":"false","fullname":"true","v-rel-max":20}
    yield "catalog", base+"cad.api?"+urlencode(cad)
    base = "https://ssd.jpl.nasa.gov/api/horizons.api?"
    for label, target, start, end in VECTOR_WINDOWS:
        params = {"format":"json","COMMAND":f"'{target};'","CENTER":"'500@599'",
                  "EPHEM_TYPE":"VECTORS","START_TIME":f"'{start}'","STOP_TIME":f"'{end}'",
                  "STEP_SIZE":"'5 d'","VEC_TABLE":2,"VEC_CORR":"NONE","REF_PLANE":"ECLIPTIC",
                  "OUT_UNITS":"KM-S","CSV_FORMAT":"YES"}
        yield label, base+urlencode(params)
    for label, target in (("oe31-sky","2019 OE31"),("uq146-sky","2025 UQ146")):
        params = {"format":"json","COMMAND":f"'{target};'","CENTER":"'500@399'",
                  "EPHEM_TYPE":"OBSERVER","START_TIME":"'2026-11-15'",
                  "STOP_TIME":"'2028-12-01'","STEP_SIZE":"'7 d'",
                  "QUANTITIES":"'1,3,9,23,36'","ANG_FORMAT":"DEG",
                  "CSV_FORMAT":"YES","EXTRA_PREC":"YES"}
        yield label, base+urlencode(params)


def table(reply):
    text = reply["result"].split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    return [[v.strip() for v in x] for x in csv.reader(text.strip().splitlines())]


def write_tables(destination, raw):
    destination.mkdir(parents=True, exist_ok=True)
    receipts = {}
    states = [["window","jd_tdb","x_km","y_km","z_km","vx_km_s","vy_km_s","vz_km_s"]]
    sky = [["target","utc","ra_deg","dec_deg","dra_cosdec_as_hr","ddec_as_hr","apmag",
            "elongation_deg","ra_3sigma_as","dec_3sigma_as"]]
    for label, url in requests():
        data = (raw/f"{label}.json").read_bytes() if raw else urlopen(url,timeout=120).read()
        receipts[label] = {"url":url,"raw_bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}
        reply = json.loads(data)
        if label == "catalog":
            (destination/"encounters.json").write_bytes(data)
        elif label.endswith("-sky"):
            for row in table(reply):
                sky.append([label.split("-")[0],*[row[i] for i in (0,3,4,5,6,7,9,11,12)]])
        else:
            for row in table(reply):
                states.append([label,row[0],*row[2:8]])
    string = io.StringIO()
    csv.writer(string,lineterminator="\n").writerows(states)
    (destination/"horizons.csv.gz").write_bytes(gzip.compress(string.getvalue().encode(),mtime=0))
    with (destination/"recovery.csv").open("w",newline="") as out:
        csv.writer(out,lineterminator="\n").writerows(sky)
    (destination/"jpl-receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--raw",type=Path,help="Use a local directory of original JSON responses")
    args = parser.parse_args()
    write_tables(args.output, args.raw)
