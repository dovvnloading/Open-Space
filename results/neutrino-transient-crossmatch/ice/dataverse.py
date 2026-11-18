"""Download the versioned Harvard releases and losslessly encode archived 90% pixel regions.

python dataverse.py [--cache PATH] [--out PATH]. The annual files are verified
against the publisher's MD5 and processed in memory; no skymap is approximated
by an enclosing box. Credentials are not required.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy_healpix import HEALPix

ROOT = Path(__file__).resolve().parent
BASE = "https://dataverse.harvard.edu/api/"
UA = "Open-Space research reproducibility (public repository github.com/dovvnloading/Open-Space)"
RELEASES = [("SCRUCD","4.0"),("RX28YT","1.0")]


def fetch(url, accept):
    req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":accept})
    for i in range(3):
        try:
            with urllib.request.urlopen(req,timeout=120) as r:
                return r.read()
        except OSError:
            if i==2:
                raise
            time.sleep(2**i)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def uniq_cells(cells, order):
    """Merge complete nested sibling quartets; preserve the pixel footprint."""
    uniq=[]
    cells=np.sort(cells)
    for level in range(order,0,-1):
        parents, counts=np.unique(cells//4,return_counts=True)
        complete=parents[counts==4]
        keep=~np.isin(cells//4,complete)
        uniq.extend((4*4**level+cells[keep]).tolist())
        cells=complete
        if not len(cells):
            break
    uniq.extend((4+cells).tolist())
    return sorted(uniq)


def legacy_fits(data, rows):
    with fits.open(io.BytesIO(gzip.decompress(data))) as f:
        h=f[1].header
        key=f'{h["RUNID"]}_{h["EVENTID"]}'
        r=rows[key]
        assert h["COORDSYS"]=="C" and h["ORDERING"]=="RING" and h["TTYPE1"]=="2LLH"
        nside=h["NSIDE"]
        order=int(np.log2(nside))
        a=f[1].data.field(0).ravel()
        assert nside==2**order and len(a)==12*nside**2
        thresholds=re.search(r"2LLH of ([\d.]+)\(([\d.]+)\)",h["COMMENTS"],re.I)
        if thresholds:
            threshold=float(thresholds[2])
        elif "based on Wilks theorem" in h["COMMENTS"]:
            # Two fitted sky coordinates: chi-square(2) 90th percentile.
            threshold=float(-2*np.log(.1))
        else:
            raise ValueError((key,h["COMMENTS"]))
        finite=np.isfinite(a)
        best=float(a[finite].min())
        selected=np.flatnonzero(finite & (a<=best+threshold))
        nested=HEALPix(nside=nside,order="ring").ring_to_nested(selected)
        cells=uniq_cells(nested,order)
        sizes=[4**(order-int(np.floor(np.log2(c/4)/2))) for c in cells]
        assert sum(sizes)==len(selected)
        return dict(run_event=key,alert_id=r["NAME"],utc=r["START"].replace(" ","T")+"Z",ra_deg=float(r["RA"]),dec_deg=float(r["DEC"]),
                    stream=r["I3TYPE"],cosmic_ray_veto=r["CR_VETO"].lower()=="true",nside=nside,delta_2llh90=threshold,
                    finite_min_2llh=best,area_deg2=float(len(selected)*4*np.pi/(12*nside**2)*(180/np.pi)**2),
                    sha256_fits_gz=sha(data),uniq90=cells)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache",type=Path,default=Path(".ice-cache"));ap.add_argument("--out",type=Path,default=ROOT)
    args=ap.parse_args()
    cache,out=args.cache,args.out
    cache.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    receipts=[]
    def obtain(name,url,accept,cache_name=None):
        p=cache/cache_name if cache_name else None
        b=p.read_bytes() if p and p.exists() else fetch(url,accept)
        if p and not p.exists():p.write_bytes(b)
        receipts.append(dict(file=name,url=url,sha256=sha(b),bytes=len(b)))
        return b
    releases=[]
    for n,(suffix,ver) in enumerate(RELEASES,1):
        url=BASE+"datasets/:persistentId/?persistentId="+urllib.parse.quote("doi:10.7910/DVN/"+suffix,safe="")
        meta=json.loads(obtain("IceCat-"+str(n)+" release metadata",url,"application/json",f"dv-meta{n}"))
        v=meta["data"]["latestVersion"]
        assert f'{v["versionNumber"]}.{v["versionMinorNumber"]}'==ver, "Release changed: review and preserve original version."
        releases.append({f["label"]:f["dataFile"] for f in v["files"]})
    def file(meta, label, cache_name=None):
        info=meta[label]
        b=obtain(label,BASE+"access/datafile/"+str(info["id"]),"application/octet-stream",cache_name)
        assert len(b)==info["filesize"]
        if not label.endswith(".tab"):
            assert hashlib.md5(b).hexdigest()==info["checksum"]["value"]
        return b
    base,rev=releases
    raw=file(base,"IceCube_Gold_Bronze_Tracks.tab","icecat1_v4.tab")
    (out/"icecat1_v4.tsv").write_bytes(raw)
    rows={r["RUNID"]+"_"+r["EVENTID"]:r for r in csv.DictReader(io.StringIO(raw.decode()),delimiter="\t")}
    assert len(rows)==348
    raw=file(rev,"Table_Preliminary_Release.tab","Table_Preliminary_Release.tab")
    (out/"icecat2_v1.tsv").write_bytes(raw)
    polygon=[]
    z=zipfile.ZipFile(io.BytesIO(file(rev,"contours.zip","contours.zip")))
    for name in sorted(z.namelist()):
        if name.startswith("contours/txt/") and name.endswith("_90.txt"):
            alert_id="IC"+re.search(r"(\d{6}[A-Z])",name)[1]
            b=z.read(name)
            xy=np.loadtxt(io.BytesIO(b),skiprows=1)
            polygon.append(dict(alert_id=alert_id,vertices_rad=xy.tolist(),sha256_contour_txt=sha(b)))
    assert len(polygon)==7
    (out/"icecat2_published90.jsonl").write_text("".join(json.dumps(p,separators=(",",":"))+"\n" for p in polygon))
    file(rev,"fits.zip","fits.zip")
    def year(label):
        tag=re.search(r"(\d{4})\.tar",label)[1]
        # Checkpoint compact masks and source receipts, never retain the large annual archives.
        checkpoint=cache/("icecat1_"+tag+"_masks.json")
        if checkpoint.exists():
            c=json.loads(checkpoint.read_text());receipts.append(c["receipt"]);return c["masks"]
        b=file(base,label,"small_2023.tar" if tag=="2023" else None)
        report=dict(file=label,url=BASE+"access/datafile/"+str(base[label]["id"]),sha256=sha(b),bytes=len(b))
        with tarfile.open(fileobj=io.BytesIO(b)) as archive:
            masks=[legacy_fits(archive.extractfile(m).read(),rows) for m in archive if m.isfile() and m.name.endswith(".fits.gz") and not m.name.rsplit("/",1)[-1].startswith("._")]
        checkpoint.write_text(json.dumps(dict(receipt=report,masks=masks),separators=(",",":")))
        return masks
    with ThreadPoolExecutor(max_workers=2) as pool:
        masks=[m for batch in pool.map(year,sorted(n for n in base if n.endswith(".tar"))) for m in batch]
    assert len(masks)==len(rows) and {r["run_event"] for r in masks}==set(rows)
    (out/"icecat1_likelihood90.jsonl").write_text("".join(json.dumps(m,sort_keys=True,separators=(",",":"))+"\n" for m in sorted(masks,key=lambda m:m["utc"])))
    with (out/"dataverse_receipts.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(receipts[0]),lineterminator="\n");w.writeheader();w.writerows(sorted(receipts,key=lambda r:r["file"]))


if __name__=="__main__":
    main()
