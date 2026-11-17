"""Compare D's catalog deficits with nearby stars on the identical eight ZTF exposures."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import zipfile

import numpy as np
from pyarrow import fs
import pyarrow.parquet as pq

import hats_fallback as hats
from d_science import DAYS
from fetch import position, separation, targets

ROOT = Path(__file__).resolve().parent
RADIUS = 120
TARGET_OID = "1544204300001642"


def api(center, exposures):
    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(dict(
        POS=f"CIRCLE {center[0]:.7f} {center[1]:.7f} {RADIUS/3600:.8f}", BANDNAME="r",
        TIME="58735 58750", FORMAT="csv", BAD_CATFLAGS_MASK=32768, COLLECTION="ztf_dr24"))
    with urlopen(url, timeout=120) as response:
        raw = response.read()
    fields = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if "expid" not in (fields.fieldnames or []):
        raise ValueError(raw[:180])
    matched = [r for r in fields if r["expid"] in exposures]
    return raw, matched, dict(url=url)


def cloud(center, exposures):
    remote = fs.S3FileSystem(anonymous=True, region="us-east-1")
    with remote.open_input_stream(hats.CAT + "/partition_info.csv") as stream:
        index = stream.read()
    leaves = {(int(r["Norder"]),int(r["Npix"])) for r in csv.DictReader(io.StringIO(index.decode()))}
    keys = [(k,p) for k in sorted({x[0] for x in leaves}) for p in hats.cone(k,*center,RADIUS) if (k,p) in leaves]
    lod=12
    cells = hats.cone(lod,*center,RADIUS)
    width=4**(29-lod)
    filters=[[("filterid","=",2),("fieldid","=",1544),("_healpix_29",">=",p*width),
              ("_healpix_29","<",(p+1)*width)] for p in cells]
    result, parts = [], []
    hmjd={k:float(v["hjd"])-2400000.5 for k,v in exposures.items()}
    for level,cell in keys:
        path = f"{hats.CAT}/dataset/Norder={level}/Dir={cell//10000*10000}/Npix={cell}.parquet"
        parts.append(dict(path="s3://"+path, bytes=remote.get_file_info(path).size))
        table = pq.read_table(path,filesystem=remote,filters=filters,columns=hats.COLUMNS)
        for obj in table.to_pylist():
            if separation(obj["objra"],obj["objdec"],*center) > RADIUS:
                continue
            series = obj["lightcurve"]
            if isinstance(series,dict):
                series=[dict(zip(series,vals)) for vals in zip(*series.values())]
            for p in series:
                for exp,time in hmjd.items():
                    if abs(p["hmjd"]-time) < 2/86400:
                        result.append(dict(expid=exp,oid=str(obj["objectid"]),ra=obj["objra"],dec=obj["objdec"],
                                           hmjd=p["hmjd"],mag=p["mag"],magerr=p["magerr"],catflags=p["catflags"]))
    stream=io.StringIO()
    w=csv.DictWriter(stream,"expid oid ra dec hmjd mag magerr catflags".split(),lineterminator="\n")
    w.writeheader(); w.writerows(result)
    return stream.getvalue().encode(),result,dict(parts=parts,index_sha256=hashlib.sha256(index).hexdigest(),radius_arcsec=RADIUS)


def main():
    t=next(t for t in targets() if t["label"]=="D")
    center=position(t,2000+(58745-51544.5)/365.25)
    with zipfile.ZipFile(ROOT / "lightcurve_inputs.zip") as archive:
        body=gzip.decompress(archive.read("cache/D.csv.gz"))
    exposures={r["expid"]:r for r in csv.DictReader(io.StringIO(body.decode())) if r["field"]=="1544" and int(float(r["mjd"])) in DAYS}
    assert len(exposures)==8
    receipt={}
    for method,fn in [("HATS",cloud),("API",api)]:
        try:
            data,rows,meta=fn(center,exposures)
            if TARGET_OID not in {str(r["oid"]) for r in rows}:
                raise ValueError("D missing in archive response")
            (ROOT/"d_catalog_control_input.csv.gz").write_bytes(gzip.compress(data,mtime=0))
            receipt[method]=dict(**meta,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),matched_points=len(rows))
            break
        except Exception as exc:
            receipt[method]=dict(error=str(exc))
            print(method,repr(exc),flush=True)
    (ROOT/"d_catalog_control_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    if all("error" in v for v in receipt.values()):
        raise ValueError("Both public catalog endpoints failed")
    by_id={}
    for r in rows:
        if int(r["catflags"])==0 and 0<float(r["magerr"])<=.1:
            by_id.setdefault(str(r["oid"]),{})[str(r["expid"])]=r
    usable={k:v for k,v in by_id.items() if set(v)==set(exposures)}
    e=list(exposures)
    other=[k for k in usable if k != TARGET_OID]
    if len(other)<3 or TARGET_OID not in usable:
        raise ValueError(f"Only {len(other)} comparison objects")
    base=[int(float(exposures[i]["mjd"]))!=58745 for i in e]
    mag=np.array([[float(usable[k][i]["mag"]) for i in e] for k in [TARGET_OID]+other])
    delta=mag-np.median(mag[:,base],axis=1)[:,None]
    med=np.median(delta[1:],axis=0)
    out=[]
    for j,exp in enumerate(e):
        spread=abs(delta[1:,j]-med[j])
        out.append(dict(expid=exp,n_controls=len(other),target_dmag=round(delta[0,j],5),
                        control_median=round(med[j],5),control_p95_abs=round(np.percentile(spread,95),5),
                        target_differential_dmag=round(delta[0,j]-med[j],5),
                        controls_fainter_than_target=sum(delta[1:,j]>=delta[0,j])))
    with (ROOT/"d_catalog_controls.csv").open("w") as stream:
        w=csv.DictWriter(stream,out[0],lineterminator="\n")
        w.writeheader(); w.writerows(out)
    print(out,flush=True)


if __name__ == "__main__":
    main()
