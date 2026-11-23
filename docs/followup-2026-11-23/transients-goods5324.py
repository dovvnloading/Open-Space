"""Bounded, receipt-preserving Stage-3 check for GOODS-S 185 in program 5324.

The adjacent protocol was frozen before examining this epoch's flux. This uses
the repository's existing measure.py aperture prescription, not a new search.
Run with --cache pointing to scratch; outputs are compact CSV/JSON receipts.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
from io import BytesIO
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROTOCOL = HERE / "transients-goods5324-protocol.json"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def run(cache, output, expected_receipts):
    protocol = json.loads(PROTOCOL.read_text())
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    expected = {}
    if expected_receipts and expected_receipts.exists():
        frozen = json.loads(expected_receipts.read_text())
        if frozen["protocol_sha256"] != digest(PROTOCOL.read_bytes()):
            raise RuntimeError("Protocol differs from the frozen receipt")
        expected = {(item["url"], r["first_byte"], r["last_byte"]): r["sha256"]
                    for item in frozen["inputs"] for r in item["ranges"]}
    receipt = {"protocol_sha256": digest(PROTOCOL.read_bytes()),
               "script_sha256": digest(Path(__file__).read_bytes()),
               "run_utc": datetime.now(timezone.utc).isoformat(),
               "base_commit": protocol["base_commit"], "inputs": []}
    coord = protocol["target"]
    ra, dec = coord["ra_deg"], coord["dec_deg"]
    gy, gx = np.mgrid[-44:45, -44:45] * .03
    rr = np.hypot(gx, gy)
    ap = rr <= .15
    ann = (rr >= .27) & (rr <= .42)
    pix_njy = np.deg2rad(.03 / 3600)**2 * 1e15
    offsets = [(x, y) for x in np.arange(-.9, .901, .3)
               for y in np.arange(-.9, .901, .3) if .45 <= np.hypot(x,y) <= 1.06]
    masks = [(np.hypot(gx-x,gy-y)<=.15,
              (np.hypot(gx-x,gy-y)>=.27)&(np.hypot(gx-x,gy-y)<=.42)) for x,y in offsets]

    def flux(image, aper, ring):
        if not np.isfinite(image[aper]).all() or np.isfinite(image[ring]).mean() < .9:
            return np.nan
        return float(np.sum(image[aper]-np.nanmedian(image[ring]))*pix_njy)

    rows = []
    for key in protocol["products"]:
        url = "https://stpubdata.s3.amazonaws.com/" + key
        band = key.split("clear-")[1].split("_")[0].upper()
        entry = {"url":url,"filter":band,"ranges":[]}
        receipt["inputs"].append(entry)
        def section(first,last):
            path=cache/(digest(f"{url}:{first}:{last}".encode())+".bin")
            cached=path.exists()
            if cached:
                raw=path.read_bytes()
            else:
                req=Request(url,headers={"Range":f"bytes={first}-{last}"})
                with urlopen(req,timeout=45) as response:
                    cr=response.headers.get("Content-Range","")
                    if response.status !=206 or not cr.startswith(f"bytes {first}-{last}/"):
                        raise RuntimeError(f"Invalid byte range: {response.status} {cr}")
                    raw=response.read()
                path.write_bytes(raw)
            if len(raw)!=last-first+1:
                raise RuntimeError("Byte-range length mismatch")
            if expected and digest(raw) != expected.get((url, first, last)):
                raise RuntimeError("Input bytes differ from the frozen receipt; use a separate campaign")
            entry["ranges"].append({"first_byte":first,"last_byte":last,"bytes":len(raw),
                                    "sha256":digest(raw),"from_cache":cached})
            return raw
        try:
            hdrs=BytesIO(section(0,57599))
            primary=fits.Header.fromfile(hdrs)
            header=fits.Header.fromfile(hdrs)
            offset=hdrs.tell()
            if (header["EXTNAME"],header["BITPIX"],header["BUNIT"]) != ("SCI",-32,"MJy/sr"):
                raise RuntimeError("Unexpected SCI format/units")
            wcs=WCS(header,naxis=2)
            nx,ny=header["NAXIS1"],header["NAXIS2"]
            x,y=wcs.all_world2pix([[ra,dec]],0)[0]
            scale=abs(header["CDELT1"])*3600
            radius=int(np.ceil(1.8/scale))
            xc,yc=int(np.rint(x)),int(np.rint(y))
            entry.update(utc=primary.get("DATE-BEG"),xy=[float(x),float(y)],shape=[ny,nx],
                         sci_data_offset=offset,crds_context=primary.get("CRDS_CTX"),
                         pipeline_version=primary.get("CAL_VER"))
            if not (radius<=xc<nx-radius and radius<=yc<ny-radius):
                entry["status"]="outside_safe_image_bounds"
                continue
            x0,y0=xc-radius,yc-radius
            y1=yc+radius+1
            raw=section(offset+y0*nx*4,offset+y1*nx*4-1)
            strip=np.frombuffer(raw,">f4").reshape(y1-y0,nx)
            stamp=strip[:,x0:xc+radius+1]
            allflux=[]
            nfinite=0
            central=None
            for dx in [-.02,0,.02]:
                for dy in [-.02,0,.02]:
                    xx,yy=wcs.all_world2pix(ra+(gx+dx)/(3600*np.cos(np.deg2rad(dec))),
                                           dec+(gy+dy)/3600,0)
                    image=map_coordinates(stamp,[yy-y0,xx-x0],order=1,mode="constant",cval=np.nan)
                    value=flux(image,ap,ann)
                    allflux.append(value)
                    if dx==dy==0:
                        blank=np.array([flux(image,a,b) for a,b in masks])
                        nfinite=int(np.isfinite(blank).sum())
                        scatter=float(1.4826*np.nanmedian(np.abs(blank-np.nanmedian(blank)))) if nfinite>=15 else np.nan
                        central=(value,scatter)
            value,scatter=central
            if not np.isfinite([value,scatter]).all() or scatter<=0:
                entry["status"]="insufficient_finite_coverage_or_noise"
                entry["finite_displaced_apertures"]=nfinite
                continue
            entry["status"]="measured_mosaic_only"
            rows.append({"site":185,"pid":5324,"obs":27,"filter":band,"utc":primary.get("DATE-BEG"),
                         "aperture_njy":value,"local_scatter_njy":scatter,
                         "offset_grid_min_njy":float(np.nanmin(allflux)),
                         "offset_grid_max_njy":float(np.nanmax(allflux)),
                         "finite_displaced_apertures":nfinite})
        except Exception as error:
            entry["status"]="failed"
            entry["error"]=str(error)
        print(band,entry["status"],flush=True)
    (output/"transients-goods5324-receipts.json").write_text(json.dumps(receipt,indent=2)+"\n")
    fields=["site","pid","obs","filter","utc","aperture_njy","local_scatter_njy",
            "offset_grid_min_njy","offset_grid_max_njy","finite_displaced_apertures"]
    # Write even an empty table, so a failed run cannot leave a stale success.
    with (output/"transients-goods5324-photometry.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    print(json.dumps(rows,indent=2))
    if any(entry["status"] != "measured_mosaic_only" for entry in receipt["inputs"]):
        raise RuntimeError("One or more required products lack a valid measurement; inspect the saved receipt")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--expected-receipts",type=Path,
                        default=HERE/"transients-goods5324-receipts.json")
    args=parser.parse_args()
    run(args.cache,args.out,args.expected_receipts)
