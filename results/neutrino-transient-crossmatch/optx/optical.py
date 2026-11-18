"""Normalize frozen survey projections. Usage: python optical.py"""
import csv
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

DIR = Path(__file__).resolve().parent
START = "2024-01-01T00:00:00Z"
END = "2026-11-17T12:09:02Z"


def iso(value):
    return datetime.fromisoformat(value.replace("Z", "")).isoformat(timespec="milliseconds") + "Z"


def jd(value, offset=0):
    utc = datetime(1970, 1, 1) + timedelta(days=float(value) + offset - 2440587.5)
    return iso(utc.isoformat())


def sex(value, hour=False):
    deg, minute, sec = map(float, value.strip().split(":"))
    sign = -1 if value.strip().startswith("-") else 1
    return sign * (abs(deg) + minute / 60 + sec / 3600) * (15 if hour else 1)


def read(name):
    return list(csv.DictReader((DIR / "data" / name).read_text().lstrip().splitlines()))


def tns_id(name):
    m = re.fullmatch(r"(?:TDE|AT|SN)?\s*((?:19|20)\d\d[a-z]+)", name.strip(), re.I)
    if m:
        return "AT" + m[1].lower()
    return name.strip()


def write(name, rows):
    with (DIR / name).open("w", newline="") as stream:
        w = csv.DictWriter(stream, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    all_rows = []
    def add(name, ra, dec, at, anchor, catalog, aliases, cls):
        if at and START[:10] <= at[:10] and at <= END:
            all_rows.append(dict(source_id=tns_id(name),ra_deg=round(float(ra),8),dec_deg=round(float(dec),8),
                                 utc=at, family="optical", anchor=anchor,catalog=catalog,aliases=aliases,classification=cls))
    for r in read("tns_selected.csv"):
        add(r["name"],r["ra"],r["declination"],iso(r["discoverydate"]),"discovery","TNS",r["internal_names"],r["type"])
    for r in read("tdecat_selected.csv"):
        add(r["AT name"],sex(r["RA"],True),sex(r["DEC"]),iso(r["Discovery date (UT)"]),"discovery","TDECat",r["ZTF name"],"TDE")
    for r in read("bts_selected.csv"):
        add(r["IAUID"] or r["ZTFID"],sex(r["RA"],True),sex(r["Dec"]),jd(r["peakt"],2458000),
            "observed_peak","ZTF BTS",r["ZTFID"] + " " + r["discID"],r["type"])
    for r in read("lions_selected.csv"):
        add(r["tns_name"] or r["ztf_name"],r["ra"],r["dec"],jd(r["peak_time_jd"]),"observed_peak",
            "Z+LIONS",r["ztf_name"]+" "+r["lsst_name"],r["classification"])
    for r in csv.DictReader((DIR / "nuclear_events.csv").open()):
        add(r["source_id"],r["ra_deg"],r["dec_deg"],r["utc"],r["anchor"],r["catalog"],"",r["classification"])
    priority = {"TNS":3,"TDECat":4,"ZTF BTS":2,"Z+LIONS":1,"Fink/ZTF Quintin":1}
    best={}
    for r in sorted(all_rows,key=lambda r:priority[r["catalog"]],reverse=True):
        best[r["source_id"]]=r
    write("optical_events.csv", sorted(best.values(),key=lambda r:(r["utc"],r["source_id"])))
    write("optical_provenance.csv", sorted(all_rows,key=lambda r:(r["source_id"],r["catalog"])))


if __name__ == "__main__":
    main()
