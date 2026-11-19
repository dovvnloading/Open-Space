#!/usr/bin/env python3
"""Join CDS COSMOS2025 photometry to a pre-frozen LePHARE row-index shortlist.

This is row-index alignment on CDS release J/A+A/704/A339, NOT a cross-release
crossmatch or astrophysical validation. The entire phot.dat source is streamed,
counted and SHA-256 hashed, and the selected bytes are parsed according to the
*same* CDS ReadMe. No whole photometry source is committed to the repo.
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import pathlib
import re
import urllib.request

BASES = (
    "https://cdsarc.cds.unistra.fr/ftp/J/A+A/704/A339/",
    "https://cdsarc.cds.unistra.fr/ftp/J/A%2BA/704/A339/",
    "https://cdsarc.u-strasbg.fr/ftp/J/A+A/704/A339/",
)
N_EXPECTED = 784016
MAX_SOURCE_BYTES = 7_000_000_000

def layout(readme, filename):
    pat = rf"Byte-by-byte Description of file:\s*{re.escape(filename)}"
    match = re.search(pat, readme, re.I)
    if not match:
        raise ValueError(f"Missing {filename} ReadMe byte layout")
    text = readme[match.end():]
    tail = re.search("Byte-by-byte Description of file:", text, re.I)
    if tail:
        text = text[:tail.start()]
    fmt = re.compile(r"^\s*(\d+)(?:-\s*(\d+))?\s+([AIFED]\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)\s*(.*)$")
    fields=[]
    for line in text.splitlines():
        m = fmt.match(line)
        if m:
            fields.append({"begin": int(m.group(1))-1, "end": int(m.group(2) or m.group(1)),
                           "format": m.group(3), "unit": m.group(4), "name": m.group(5),
                           "description": m.group(6)})
        elif fields and line.strip():
            fields[-1]["description"] += " "+line.strip()
    assert len(fields) >= 200, "Unrecognized phot.dat layout"
    assert max(x["end"] for x in fields) >= 5500, "Truncated phot.dat layout"
    return fields

def download_readme():
    errors=[]
    for base in BASES:
        try:
            with urllib.request.urlopen(urllib.request.Request(base+"ReadMe",headers={"User-Agent":"Open-Space/1.0 (research)"}), timeout=120) as r:
                raw=r.read()
            if len(raw)<1000:
                raise ValueError("too few ReadMe bytes")
            return raw, base+"ReadMe"
        except Exception as exc:
            errors.append(str(exc))
    raise RuntimeError("No official CDS ReadMe: "+" | ".join(errors))

def get_gzip_stream():
    errors=[]
    for base in BASES:
        for name in ("phot.dat.gz","phot.dat"):
            url=base+name
            try:
                response=urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"Open-Space/1.0 (research)"}),timeout=120)
                return response,url,name.endswith(".gz")
            except Exception as exc:
                errors.append((url,str(exc)))
    raise RuntimeError("Cannot access full CDS photometry file: "+str(errors))

def selected_fields(layout):
    base = {"Id","Tile","Segment-Id","Id-spz-K25","RAdeg","DEdeg","Ximage","Yimage"}
    # Keep all supplied 37-band model fluxes AND their error/quality columns,
    # while excluding large, unrelated morphology/aper subdivisions.
    chosen=[x for x in layout if x["name"] in base
            or x["name"].startswith(("Flux-mod-","e_Flux-mod-","e_Flux-u-mod-"))
            or x["name"].lower() in ("flag", "flags", "flag-miri", "flag-nircam")]
    actual={x["name"] for x in chosen}
    if not {"Id","RAdeg","DEdeg","Tile"}.issubset(actual):
        raise ValueError("Required source association columns are unavailable")
    if len(chosen)<35:
        raise ValueError("Fewer than 35 photometric/association columns available")
    return chosen

def parse(line, f):
    raw=line[f["begin"]:f["end"]].strip()
    if not raw:
        return None
    if f["format"].startswith("A"):
        return raw
    try:
        val=float(raw.replace("D","E"))
    except ValueError:
        return None
    import math
    if not math.isfinite(val) or val in (-9999.,-999.9,-999.,-99.9,-99.,999.,9999.):
        return None
    return val

def run(out, shortlist):
    with shortlist.open(newline="") as f:
        targets={int(r["source_row_1based"]): r for r in csv.DictReader(f)}
    if not targets:
        raise ValueError("Shortlist CSV contains no rows")
    if any(not 1<=i<=N_EXPECTED for i in targets):
        raise ValueError("Shortlist contains invalid row indices")
    raw,readme_url=download_readme()
    selected=selected_fields(layout(raw.decode("latin-1"),"phot.dat"))
    started=dt.datetime.now(dt.timezone.utc).isoformat()
    response,url,compressed=get_gzip_stream()
    sha=hashlib.sha256()
    count=0
    total=0
    matched={}
    try:
        # Hash exactly the REMOTE file bytes while streaming the decompressed
        # 5,563-byte records. Limit decompressed bytes and validate rows.
        import io
        class HashStream(io.RawIOBase):
            def __init__(self,r):
                self.r=r
            def readable(self):
                return True
            def readinto(self,b):
                nonlocal total
                block=self.r.read(len(b))
                if block:
                    sha.update(block)
                    total+=len(block)
                    if total>MAX_SOURCE_BYTES:
                        raise ValueError("remote input exceeds stream safety bound")
                    b[:len(block)]=block
                return len(block)

        # gzip.GzipFile accepts fileobj with a read() method; decompressed
        # data iterator is exposed via io.TextIOWrapper.
        wrapped=io.BufferedReader(HashStream(response),buffer_size=1024*1024)
        source=gzip.GzipFile(fileobj=wrapped) if compressed else wrapped
        with io.TextIOWrapper(source,encoding="ascii",errors="strict") as inp:
            for line in inp:
                if not line.strip():
                    continue
                count+=1
                if count in targets:
                    record={x["name"]:parse(line,x) for x in selected}
                    record["source_row_1based"]=count
                    record["zfinal_marginal"]=targets[count].get("zfinal")
                    record["mass_med_marginal"]=targets[count].get("mass_med")
                    record["mass_u68_marginal"]=targets[count].get("mass_u68")
                    matched[count]=record
    finally:
        response.close()
    if count!=N_EXPECTED:
        raise ValueError(f"Expected {N_EXPECTED} photometry rows; got {count}; do not interpret partial source")
    if set(matched)!=set(targets):
        raise ValueError("Lost selected row(s) during photometry join")
    ids=[v["Id"] for v in matched.values()]
    if len(ids)!=len(set(ids)):
        raise ValueError("Photometry Id not unique within candidate set")
    for v in matched.values():
        if not (140 < v["RAdeg"] < 160 and 0 < v["DEdeg"] < 4):
            raise ValueError("A COSMOS source has implausible J2000 coordinates")
    out.mkdir(parents=True,exist_ok=True)
    with (out/"highz_candidate_photometry.csv").open("w",newline="") as f:
        names=["source_row_1based"]+[x["name"] for x in selected]+[
            "zfinal_marginal","mass_med_marginal","mass_u68_marginal"]
        writer=csv.DictWriter(f,fieldnames=names)
        writer.writeheader()
        for _,row in sorted(matched.items()):
            writer.writerow(row)
    receipt={"started_utc":started,"completed_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
             "photometry_url":url,"file_compressed":compressed,"sha256_raw_download":sha.hexdigest(),
             "raw_source_bytes":total,"readme_url":readme_url,
             "readme_sha256":hashlib.sha256(raw).hexdigest(),"record_count":count,
             "joined_shortlist_rows":len(matched),"columns":list(matched[next(iter(matched))]),
             "release_warning":"CDS mirror not demonstrated equivalent to official 2026-May-04 IAP v1.1.",
             "posterior_warning":"Photometry is NOT a second SED refit or a joint p(M*,z) posterior."}
    (out/"photometry_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))
    return receipt

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shortlist",type=pathlib.Path,required=True)
    p.add_argument("--out",type=pathlib.Path,required=True)
    a=p.parse_args()
    run(a.out,a.shortlist)

if __name__=="__main__":
    main()
