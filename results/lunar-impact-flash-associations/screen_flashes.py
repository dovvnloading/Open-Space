"""Replay the spatial screen against the NASA flash and LROC impact catalogues."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent
URLS = {
    "craters.zip": "https://pds.lroc.im-ldi.com/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/EXTRAS/SHAPEFILE/GLOBAL_NEW_CRATERS/GLOBAL_NEW_CRATERS_180.ZIP",
    "flash.pdf": "https://www.nasa.gov/wp-content/uploads/2016/12/alamo_lunar_impact_observations443.pdf",
    "pos.pdf": "https://ntrs.nasa.gov/api/citations/20150021386/downloads/20150021386.pdf",
}

def distance(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(math.radians, (lon1, lat1, lon2, lat2))
    h = math.sin((lat2-lat1)/2)**2 + math.cos(lat1)*math.cos(lat2)*math.sin((lon2-lon1)/2)**2
    return 2*1737.4*math.asin(min(1, math.sqrt(h)))

def load_craters(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        dbf = z.read("GLOBAL_NEW_CRATERS_180.DBF")
    count, header, size = struct.unpack_from("<IHH", dbf, 4)
    fields = [(dbf[i:i+11].rstrip(b"\0").decode(), dbf[i+16]) for i in range(32, header-1, 32)]
    rows = []
    for index in range(count):
        pos = header + index*size + 1
        row = {}
        for key, width in fields:
            row[key] = dbf[pos:pos+width].decode().strip()
            pos += width
        rows.append(row)
    return rows

def pdf_text(path):
    return subprocess.check_output(["pdftotext", "-layout", str(path), "-"]).decode().replace("–", "-")

def load_flashes(cache):
    flashes = {}
    pat = re.compile(r"^\s*(\d+)\s+(\d+-\w{3}-\d{4})\s+([\d:.]+).*?([\d.]+)\s*([EW])\s+([\d.]+)\s*([NS])\b")
    for line in pdf_text(cache/"flash.pdf").splitlines():
        m = pat.search(line)
        if m:
            idx, date, time, lon, ew, lat, ns = m.groups()
            flashes[int(idx)] = (dt.datetime.strptime(date, "%d-%b-%Y").date(), float(lon)*(-1 if ew=="W" else 1), float(lat)*(-1 if ns=="S" else 1))
    for line in pdf_text(cache/"pos.pdf").splitlines():
        cols = line.split()
        k = next((i for i,s in enumerate(cols[:3]) if re.fullmatch(r"\d{4}-\w{3}-\d{2}",s)),None)
        if k and cols[0].isdigit() and len(cols)>=k+11 and cols[k+6]!="***":
            flashes[int(cols[0])] = (dt.datetime.strptime(cols[k],"%Y-%b-%d").date(),float(cols[k+6]),float(cols[k+7]))
    return flashes

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    args=parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    manifest = json.loads((ROOT/"source_hashes.json").read_text())
    for name,url in URLS.items():
        path = args.cache/name
        if not path.exists():
            with urllib.request.urlopen(url,timeout=120) as response:
                path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest[name]:
            raise ValueError(f"Source changed: {name}")
    craters=load_craters((args.cache/"craters.zip").read_bytes())
    flashes=load_flashes(args.cache)
    dates={r["product"]:r["date"] for r in csv.DictReader((ROOT/"image_dates.csv").open())}
    print(f"{len(craters)} crater entries, {len(flashes)} positioned flashes")
    for row in craters:
        for idx,(date,lon,lat) in flashes.items():
            km=distance(float(row["lon"]),float(row["lat"]),lon,lat)
            if km < 25:
                first,last=[dates[row[k]] for k in ("before_img","after_img")]
                print(idx, f"{km:.2f} km", date, first, last, first<=date.isoformat()<=last)
    date,lon,lat=flashes[296]
    print("Independent 2013-03-17 control at 20.714 N, 24.330 W:", round(distance(-24.330,20.714,lon,lat),2), "km")
if __name__=="__main__":
    main()
