"""Acquire a separately resolved pre-ZTF PTF survey light curve, with a hard audit trail.

No success is inferred from endpoint availability: images A and B must each
have independently calibrated, clean measurements at their fixed Gaia positions.
Output is saved as a downloadable CI artifact for offline scrutiny.
"""
import csv, hashlib, io, json, sys, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "runner-output"
OUT.mkdir(parents=True, exist_ok=True)
RA, DEC, RADIUS = 126.5341, 70.0456, 0.004
POS = {"A": (126.53489916666669,70.04488027777778),
       "B": (126.53333583333334,70.04640916666666)}

def request(url):
    req = urllib.request.Request(url, headers={"User-Agent":"Open-Space lens extension 2026 contact by GitHub repo"})
    with urllib.request.urlopen(req, timeout=55) as stream:
        return stream.read()

def parse_votable(raw):
    root = ET.fromstring(raw)
    fields = [x.attrib.get("name","") for x in root.findall(".//{*}FIELD")]
    rows = [[(cell.text or "").strip() for cell in tr.findall("./{*}TD")]
            for tr in root.findall(".//{*}TR")]
    if not fields:
        raise ValueError("No VOTable FIELD declarations")
    if any(len(row) != len(fields) for row in rows):
        raise ValueError("Inconsistent VOTable row length")
    return fields,rows

def cone(table, fmt="votable"):
    query = urllib.parse.urlencode({"table":table,"RA":RA,"DEC":DEC,"SR":RADIUS,"format":fmt})
    url = "https://irsa.ipac.caltech.edu/SCS?"+query
    raw = request(url)
    (OUT / (table+"."+("xml" if fmt=="votable" else "txt"))).write_bytes(raw)
    print("SOURCE",url, "BYTES",len(raw), "SHA256",hashlib.sha256(raw).hexdigest())
    if fmt=="votable":
        return parse_votable(raw)
    return [],[]

def ptf_extract():
    for table in ("ptf_objects","ptf_lightcurves"):
        try:
            fields,rows=cone(table)
            print("PTF",table,"N",len(rows),"COLS",fields[:18], "ERRORS", " ".join((OUT/(table+".xml")).read_text(errors="replace").split("INFO")[-2:])[:260] if not rows else "")
            if not rows: continue
            d=[dict(zip(fields,row)) for row in rows]
            seen=Counter((v.get("oid"),v.get("fid")) for v in d)
            print("PTF_IDS",seen.most_common(20))
            if table=="ptf_lightcurves":
                keys=["oid","ra","dec","fid","obsmjd","mag_autocorr","magerr_auto","goodflag","photcalflag","flags","qaflags","limitmag"]
                with (OUT/"ptf-selected.csv").open("w",newline="") as f:
                    w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
                    for row in d:w.writerow({k:row.get(k,"") for k in keys})
                print("PTF_SAMPLE",json.dumps([{k:r.get(k) for k in keys} for r in d[:5]]))
        except Exception as exc:
            print("PTF_ERROR",table,type(exc).__name__,str(exc)[:450])
    # Alternative TAP query, selecting just the columns required; avoids too-wide VOTable.
    adql=("SELECT TOP 12000 oid,ra,dec,fid,obsmjd,mag_autocorr,magerr_auto,goodflag,"
          "photcalflag,flags,qaflags,limitmag FROM ptf_lightcurves WHERE "
          "1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',126.5341,70.0456,0.004))")
    url="https://irsa.ipac.caltech.edu/TAP/sync?"+urllib.parse.urlencode(
        {"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":adql})
    try:
        raw=request(url)
        (OUT/"ptf-tap.csv").write_bytes(raw)
        rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        print("TAP_BYTES",len(raw),"SHA256",hashlib.sha256(raw).hexdigest(),
              "ROWS",len(rows),"HEADER",list(rows[0]) if rows else raw[:500].decode(errors="replace"))
        print("TAP_IDS",Counter((r.get("oid"),r.get("fid")) for r in rows).most_common(20))
        print("TAP_SAMPLE",json.dumps(rows[:2])[:1500])
    except Exception as exc:
        print("TAP_ERROR",type(exc).__name__,str(exc)[:600])

def export_existing():
    root=HERE.parent
    for name in ("ztf-nights.csv","gaia-nights.csv","ps1-detections.csv","targets.csv"):
        data=(root/name).read_bytes()
        if name=="ztf-nights.csv":
            lines=data.decode().splitlines()
            data=("\n".join([lines[0]]+[l for l in lines[1:] if l.startswith("J0826+7002,")])+"\n").encode()
        (OUT/name).write_bytes(data)
        print("EXPORTED",name,len(data),hashlib.sha256(data).hexdigest())

if __name__=="__main__":
    export_existing()
    ptf_extract()
