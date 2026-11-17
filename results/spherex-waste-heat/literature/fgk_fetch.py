"""Fetch the CatWISE and NEOWISE neighborhoods of all 53 FGK candidates."""
import argparse
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
TABLES = {
    "catwise": ("catwise_2020", "source_name,ra,dec,ra_pm,dec_pm,meanobsmjd,ab_flags,cc_flags,w1mpro_pm,w1sigmpro_pm,w2mpro_pm,w2sigmpro_pm,w1snr_pm,w2snr_pm", 2015.405),
    "neo": ("neowiser_p1bs_psd", "source_id,scan_id,frame_num,ra,dec,mjd,qual_frame,qi_fact,saa_sep,cc_flags,moon_masked,ph_qual,nb,w1sat,w2sat,w1mpro,w1sigmpro,w1snr,w1rchi2,w2mpro,w2sigmpro,w2snr,w2rchi2", 2024.25),
}

def position(row, epoch):
    elapsed = epoch - float(row["ref_epoch"])
    dec = float(row["dec"])
    return (float(row["ra"]) + elapsed * float(row["pmra"]) / (3.6e6 * math.cos(math.radians(dec))),
            dec + elapsed * float(row["pmdec"]) / 3.6e6)

def read_result(kind, data):
    if kind == "catwise":
        reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
        if reader.fieldnames != TABLES[kind][1].split(","):
            raise ValueError(data[:600])
        return list(reader)
    root = ET.fromstring(data)
    bad = [el for el in root.iter() if el.tag.split("}")[-1] == "INFO" and el.attrib.get("value") in ("ERROR","OVERFLOW")]
    if bad:
        raise ValueError(ET.tostring(bad[0]))
    tab = root.find(".//{*}TABLE")
    if tab is None:
        raise ValueError(data[:400])
    columns = [el.attrib["name"] for el in tab.findall("{*}FIELD")]
    if not set(TABLES[kind][1].split(",")).issubset(columns):
        raise ValueError(f"Unexpected VOTable columns {columns}")
    return [dict(zip(columns, [el.text or "" for el in row.findall("{*}TD")], strict=True))
            for row in tab.findall(".//{*}TR")]

def fetch(target, kind, cache):
    sid = target["source_id"]
    table, columns, epoch = TABLES[kind]
    ra, dec = position(target, epoch)
    query = f"SELECT {columns} FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.0033333333))=1"
    if kind == "neo":
        query += " AND mjd>=60291 AND mjd<60524"
        url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?" + urlencode(dict(
            catalog=table, spatial="Cone", objstr=f"{ra:.8f} {dec:.8f}", radius=12, radunits="arcsec",
            selcols=columns, constraints="mjd>=60291 AND mjd<60524", outfmt=3, outrows=20000))
    else:
        url = TAP + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=20000, QUERY=query))
    base = cache / f"{kind}-{sid}"
    body = base.with_suffix(".xml" if kind == "neo" else ".csv")
    receipt = base.with_suffix(".json")
    if receipt.exists() and body.exists():
        meta = json.loads(receipt.read_text())
        if meta["url"] == url and meta["sha256"] == hashlib.sha256(body.read_bytes()).hexdigest():
            return
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={"User-Agent":"Open-Space astronomical research"}), timeout=90) as response:
                data = response.read()
            n = len(read_result(kind, data))
            if n >= 20000:
                raise ValueError("Archive row limit reached")
            meta = dict(source_id=sid, query=query, url=url, rows=n, bytes=len(data),
                        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), sha256=hashlib.sha256(data).hexdigest())
            body.write_bytes(data)
            receipt.write_text(json.dumps(meta))
            print(kind,sid,n,flush=True)
            return
        except Exception as exc:
            if attempt == 2:
                raise RuntimeError(f"{kind} {sid}") from exc
            time.sleep(2**attempt)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cache",type=Path,default=Path("fgk-cache"))
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--kind",choices=TABLES,required=True)
    a=ap.parse_args()
    a.cache.mkdir(exist_ok=True,parents=True)
    with (ROOT/"fgk_astrometry.csv").open() as stream:
        rows=list(csv.DictReader(stream))
    with concurrent.futures.ThreadPoolExecutor(a.workers) as pool:
        jobs=[pool.submit(fetch,row,a.kind,a.cache) for row in rows]
        for job in concurrent.futures.as_completed(jobs):
            job.result()
    records, receipts = [], []
    for row in rows:
        name=row["source_id"]
        base=a.cache / f"{a.kind}-{name}"
        data=base.with_suffix(".xml" if a.kind == "neo" else ".csv").read_bytes()
        records.append(json.dumps({"source_id":name,"response":data.decode("utf-8")}))
        receipts.append(json.loads(base.with_suffix(".json").read_text()))
    (ROOT/f"fgk_{a.kind}_raw.jsonl.gz").write_bytes(gzip.compress(("\n".join(records)+"\n").encode(),mtime=0))
    (ROOT/f"fgk_{a.kind}_receipts.json").write_text(json.dumps(receipts,indent=1)+"\n")

if __name__ == "__main__":
    main()
