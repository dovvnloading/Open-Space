"""Look for public Hubble images covering the complete unlisted configuration."""
import csv
import hashlib
import io
import json
import re
import sys
import time
import zipfile
from xml.etree import ElementTree
from pathlib import Path

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from PIL import Image


def fetch(url, params, receipts, dest=None):
    for attempt in range(2):
        try:
            r=requests.get(url,params=params,timeout=(12,55))
            receipts.append(dict(url=r.url,http=r.status_code,bytes=len(r.content),
                                 sha256=hashlib.sha256(r.content).hexdigest()))
            r.raise_for_status()
            if dest:
                dest.write_bytes(r.content)
            return r
        except (requests.RequestException, OSError) as exc:
            if attempt or (getattr(exc,'response',None) is not None and exc.response.status_code in (400,404,429)):
                raise
            time.sleep(2)


def covers(stcs,point):
    """Test simple STC-S circles/polygons and their unions. Unknown stays unknown."""
    s=(stcs or "").upper()
    if not s or re.search(r"\b(NOT|INTERSECTION|DIFFERENCE|BOX|ELLIPSE)\b",s):
        return None
    parts=re.split(r'\b(POLYGON|CIRCLE)\b',s)[1:]
    if not parts:
        return None
    valid=False
    for shape,part in zip(parts[::2],parts[1::2]):
        v=[float(a) for a in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:E[-+]?\d+)?',part)]
        if shape=="CIRCLE" and len(v)==3:
            valid=True
            if point.separation(SkyCoord(v[0],v[1],unit="deg")).deg<=v[2]:
                return True
        elif shape=="POLYGON" and len(v)>=6 and len(v)%2==0:
            # Straight edges of a convex spherical detector polygon remain
            # straight in the gnomonic projection, also across the RA meridian.
            lon,lat=point.spherical_offsets_to(SkyCoord(v[::2],v[1::2],unit="deg"))
            if np.any(np.cos(lon.radian)*np.cos(lat.radian)<=0):
                continue
            x=np.tan(lon.radian)
            y=np.tan(lat.radian)/np.cos(lon.radian)
            k=np.flatnonzero((y>0)!=(np.roll(y,1)>0))
            j=(k-1)%len(y)
            crossing=x[k]-y[k]*(x[j]-x[k])/(y[j]-y[k])
            valid=True
            if np.count_nonzero(crossing>0)%2:
                return True
    return False if valid else None


def caom(coords,history):
    # MAST retrieval is owned by object_history.py. Its cone matches alone do
    # not prove that either the inferred mass centre or both images are exposed.
    rows=[]
    for x in history.get("high_resolution_observations",[]):
        if x.get("obs_collection") not in ("HST","HLA") or x.get("dataproduct_type")!="image":
            continue
        coverage=[covers(x.get("s_region"),p) for p in coords]
        rows.append(dict(obs_id=x.get("obs_id"),instrument=x.get("instrument_name"),
                         filters=x.get("filters"),proposal=x.get("proposal_id"),
                         product=x.get("dataURL"),s_region=x.get("s_region"),coverage=coverage))
    full=[x for x in rows if all(v is True for v in x["coverage"])]
    unknown=sum(any(v is None for v in x["coverage"]) for x in rows)
    status=("covers_centre_and_all_images" if full else
            "unavailable_or_incomplete" if history.get("status")!="complete" or unknown else
            "no_covering_HST_image_in_completed_MAST_query")
    return dict(status=status,image_records=len(rows),unknown_footprints=unknown,covering=full,
                query_receipt=history.get("receipt"),footprints=rows)


def sources(coords,folder,receipts):
    coord=coords[0]
    r=fetch("https://catalogs.mast.stsci.edu/api/v0.1/hsc/v3/summary/magauto.csv",
            dict(ra=coord.ra.deg,dec=coord.dec.deg,radius=6/3600,pagesize=5000),receipts,folder/"hsc.csv")
    reader=csv.DictReader(io.StringIO(r.content.decode("utf-8-sig")))
    names={k.lower():k for k in reader.fieldnames or []}
    if not {"matchra","matchdec","matchid"}<=names.keys():
        raise ValueError("HSC catalog columns not present")
    rows=[]
    for entry in reader:
        p=SkyCoord(float(entry[names["matchra"]]),float(entry[names["matchdec"]]),unit="deg")
        if coords[0].separation(p).arcsec<=6:
            rows.append(dict(matchid=entry[names["matchid"]],ra=p.ra.deg,dec=p.dec.deg,
                             separations_to_centre_and_images_arcsec=coords.separation(p).arcsec.tolist()))
    return dict(status="query_complete",sources=rows)


def retain(payload,coords,prefix):
    with fits.open(io.BytesIO(payload),memmap=False) as hdus:
        for hdu in hdus:
            if hdu.data is None or hdu.data.ndim!=2 or hdu.name.upper() not in ("PRIMARY","SCI","CUTOUT"):
                continue
            try:
                w=WCS(hdu.header,hdus).celestial
                sc=float(np.mean(proj_plane_pixel_scales(w))*3600)
                if not .01<sc<.20:
                    continue
                field=Cutout2D(np.asarray(hdu.data,dtype='f4'),coords[0],(12*u.arcsec,12*u.arcsec),wcs=w,
                               mode="partial",fill_value=np.nan,copy=True)
                pix=np.column_stack(field.wcs.world_to_pixel(coords))
                if not np.all(np.isfinite(pix)):
                    continue
                yy,xx=np.indices(field.data.shape)
                for px,py in pix:
                    core=field.data[(xx-px)**2+(yy-py)**2<(.28/sc)**2]
                    if len(core)<6 or np.mean(np.isfinite(core)&(core!=0))<.8 or np.nanstd(core)==0:
                        break
                else:
                    arr=field.data.astype('f4')
                    hdr=field.wcs.to_header(relax=True)
                    for key in ('BUNIT','FILTER','FILTER1','FILTER2','INSTRUME','EXPTIME','PHOTFLAM','PHOTPLAM'):
                        if key in hdu.header or key in hdus[0].header:
                            hdr[key]=hdu.header.get(key,hdus[0].header.get(key))
                    fits.PrimaryHDU(arr,header=hdr).writeto(prefix.with_suffix('.fits'),overwrite=True)
                    med=np.nanmedian(arr)
                    noise=max(float(np.nanmedian(abs(arr-med))*1.4826),1e-9)
                    scaled=np.arcsinh((arr-med)/noise)
                    low,high=np.nanpercentile(scaled,[1,99.5])
                    gray=(255*np.clip(np.nan_to_num((scaled-low)/max(high-low,.01)),0,1)).astype('uint8')
                    Image.fromarray(gray[::-1]).save(prefix.with_suffix('.png'))
                    return dict(fits=prefix.with_suffix('.fits').name,png=prefix.with_suffix('.png').name,
                                pixel_scale_arcsec=sc,centre_and_image_xy=pix.tolist(),
                                sha256=hashlib.sha256(prefix.with_suffix('.fits').read_bytes()).hexdigest())
            except (TypeError,ValueError,IndexError):
                continue
    raise ValueError("FITS lacks confirmed usable pixels at centre and all components")


def inspect(packet,folder,history):
    folder.mkdir(exist_ok=True,parents=True)
    centre=packet['imaging']['lens_model']['center_radec_deg']
    components=packet['components']
    coords=SkyCoord([centre[0]]+[c['ra'] for c in components],
                    [centre[1]]+[c['dec'] for c in components],unit="deg")
    ra,dec=centre
    receipts=[]
    result=dict(id=packet["system_id"],centre=centre,verified_stamps=[],failures=[],inputs=receipts)
    try:
        result["caom"]=caom(coords,history)
    except Exception as exc:
        result["caom"]=dict(status="unavailable",error=str(exc))
    try:
        result["hsc"]=sources(coords,folder,receipts)
    except Exception as exc:
        result["hsc"]=dict(status="unavailable",error=str(exc))
    # A positional SIA result may merely overlap the search cone. Never call
    # its source covered until actual non-blank FITS pixels pass the WCS test.
    datasets=[]
    try:
        r=fetch("https://hla.stsci.edu/cgi-bin/hlaSIAP.cgi",
                dict(POS=f"{ra},{dec}",SIZE=.002,FORMAT="image/fits"),receipts,folder/"hla.vot")
        root=ElementTree.fromstring(r.content)
        for elt in root.iter():
            if (elt.attrib.get("name","").upper()=="QUERY_STATUS"
                    and elt.attrib.get("value","").upper() in ("ERROR","OVERFLOW")):
                raise ValueError("HLA SIA query "+elt.attrib["value"]+": "+(elt.text or "")[:130])
        table=Table.read(io.BytesIO(r.content),format="votable")
        table.convert_bytestring_to_unicode()
        result["hla"]=dict(status="query_complete",regional_rows=len(table))
        pattern=r"\bhst_\d{4,6}_[a-z0-9]+_(?:acs|wfc3|wfpc2|nicmos)_[a-z0-9_]+"
        for row in table:
            for name in table.colnames:
                for match in re.findall(pattern,str(row[name]).lower()):
                    datasets.append(re.sub(r'_(?:drz|drc|sci)$','',match))
    except Exception as exc:
        result["hla"]=dict(status="unavailable",error=str(exc))
    caom_rows=result.get("caom",{}).get("covering",[])
    datasets+= [x["obs_id"] for x in caom_rows if str(x.get("obs_id","")).lower().startswith("hst_")]
    for name in list(dict.fromkeys(datasets))[:12]:
        if len(result["verified_stamps"])>=6:
            break
        try:
            r=fetch("https://hla.stsci.edu/cgi-bin/fitscut.cgi",
                    dict(ra=ra,dec=dec,red=name,size=300,format="fits"),receipts)
            stamp=retain(r.content,coords,folder/("hla-"+str(len(result["verified_stamps"]))))
            result["verified_stamps"].append(dict(dataset=name,**stamp))
        except Exception as exc:
            result["failures"].append(dict(dataset=name,error=str(exc)))
    if not result["verified_stamps"]:
        # A covering individual/public drizzled product can be available even
        # when HLA has no matching combined dataset or HAP mosaic.
        for source in caom_rows[:6]:
            uri=str(source.get("product") or "")
            if not uri.lower().startswith("mast:hst/product/") or not re.search(
                    r"_(?:drz|drc|drw|c0m|flt|flc)\.fits(?:\.gz)?$",uri,re.I):
                continue
            try:
                r=fetch("https://mast.stsci.edu/api/v0.1/Download/file",dict(uri=uri),receipts)
                stamp=retain(r.content,coords,folder/("hst-"+str(len(result["verified_stamps"]))))
                result["verified_stamps"].append(dict(dataset=uri,**stamp))
            except Exception as exc:
                result["failures"].append(dict(dataset=uri,error=str(exc)))
            if len(result["verified_stamps"])>=3:
                break
    if not result["verified_stamps"]:
        try:
            r=fetch("https://mast.stsci.edu/hapcut/api/v0.1/astrocut",
                    dict(ra=ra,dec=dec,x=12/3600,y=12/3600,units="d"),receipts)
            with zipfile.ZipFile(io.BytesIO(r.content)) as archive:
                for member in archive.infolist():
                    if not member.filename.lower().endswith(".fits") or member.file_size>15000000:
                        continue
                    try:
                        stamp=retain(archive.read(member),coords,folder/("hap-"+str(len(result["verified_stamps"]))))
                        result["verified_stamps"].append(dict(dataset=member.filename,**stamp))
                    except ValueError as exc:
                        result["failures"].append(dict(dataset=member.filename,error=str(exc)))
                    if len(result["verified_stamps"])>=6:
                        break
        except Exception as exc:
            result["failures"].append(dict(service="HAPcut",error=str(exc)))
    return result


def run(evidence,out,history_path=None):
    evidence,out=Path(evidence),Path(out)
    out.mkdir(parents=True,exist_ok=True)
    unique={}
    for line in (evidence/'review_packets.jsonl').read_text().splitlines():
        p=json.loads(line)
        if not p['prior'] and p['imaging']['imaging_model_pass']:
            key=tuple(sorted(str(c["source_id"]) for c in p['components']))
            unique.setdefault(key,p)
    history={}
    if history_path and Path(history_path).exists():
        obj=json.loads(Path(history_path).read_text())
        history={v['system_id']:v['mast'] for v in obj['objects']}
    rows=[inspect(p,out/p["system_id"],history.get(p["system_id"],{})) for p in unique.values()]
    (out/'coverage.json').write_text(json.dumps(dict(systems=rows,scope="MAST public HST and HLA image "
        "footprints near the field and checksummed HLA/HAP cutouts. HSC alone is not an exposure inventory. "
        "An unavailable endpoint, unparsed footprint, or unusable cutout is not a non-detection."),indent=2)+"\n")


if __name__=='__main__':
    run(*sys.argv[1:])
