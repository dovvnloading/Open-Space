"""Freeze the release inputs, inspect their schemas and identify prior surveys."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import os

import requests
from astropy.io import fits
from astropy.table import Table


DOC = Path("results/xmmsl2-erass1-rass/release")
DOC.mkdir(parents=True, exist_ok=True)
CACHE = Path("work/downloads")
CACHE.mkdir(parents=True, exist_ok=True)
HEA = "https://heasarc.gsfc.nasa.gov/FTP/"
ERO = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/"
URLS = {
    "xmm_clean.fits.gz": HEA+"xmm/data/catalogues/xmmsl2_clean.fits.gz",
    "xmm_full.fits.gz": HEA+"xmm/data/catalogues/xmmsl2_total.fits.gz",
    "ero_main.tar.gz": ERO+"eRASS1_Main.v1.2.fits.tar.gz",
    "ero_supp.tar.gz": ERO+"eRASS1_Supp.tar.gz",
    "rass.fits.gz": HEA+"heasarc/dbase/tdat_files/heasarc_rass2rxs.tdat.gz",
    "pspc.rsp":HEA+"caldb/data/rosat/pspc/cpf/matrices/pspcb_gain1_256.rsp",
}
# The HEASARC FITS export is delivered by the VO TAP service.
URLS["rass.fits.gz"] = HEA+"heasarc/dbase/msc/heasarc_rass2rxs.fits.gz"


def inspect_source(item):
    key, url = item
    target = CACHE/key
    rec = {"file": key, "url":url}
    try:
        target.unlink(missing_ok=False) if target.exists() and target.stat().st_size==0 else None
        if not target.exists():
            with requests.get(url, stream=True, timeout=(12,150)) as response:
                response.raise_for_status()
                with target.with_suffix(".part").open("wb") as fp:
                    for chunk in response.iter_content(1024*1024):
                        fp.write(chunk)
                target.with_suffix(".part").replace(target)
        with target.open("rb") as fp:
            rec["sha256"] = hashlib.file_digest(fp,"sha256").hexdigest()
        rec["bytes"] = target.stat().st_size
        with tempfile.TemporaryDirectory() as tmp:
            path = target
            if key.endswith(".tar.gz"):
                with tarfile.open(target,"r|gz") as tf:
                    for member in tf:
                        if member.isfile() and member.name.lower().endswith(".fits"):
                            from shutil import copyfileobj
                            path = Path(tmp)/"table.fits"
                            with path.open("wb") as stream:
                                copyfileobj(tf.extractfile(member),stream)
                            break
            with fits.open(path,memmap=False if path.suffix==".gz" else True,lazy_load_hdus=True) as hdus:
                hdu=hdus[1]
                rec["rows"] = hdu.header.get("NAXIS2")
                rec["columns"] = [(str(c.name), str(c.unit or ""),str(c.format)) for c in hdu.columns]
                visible=[c.name for c in hdu.columns if c.format not in ("PE","PI","PJ")][:180]
                field = hdu.data[0]
                rec["example"] = {k:str(field[k])[:120] for k in visible}
    except Exception as exc:
        rec["error"] = repr(exc)
    (DOC/(key+".json")).write_text(json.dumps(rec,indent=2)+"\n")
    print(key,rec.get("rows"),rec.get("error",""),flush=True)


def metadata():
    queries={
        "rass_field":"https://heasarc.gsfc.nasa.gov/FTP/rosat/data/pspc/processed_data/rass/release/rs930101n00/",
        "s3_rass":"https://nasa-heasarc.s3.amazonaws.com/?list-type=2&prefix=rosat/data/pspc/processed_data/rass/release/&max-keys=36",
        "tns_http":"https://www.wis-tns.org/search?format=csv&ra=244.0009&decl=22.268&radius=20&coords_unit=arcsec",
        "rassmaster":"https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync",
        "eul_live":"https://erosita.mpe.mpg.de/erodat/upperlimit/service?ra=194.0&dec=-5.7&band=024&dr_survey=DR1_eRASS1",
    }
    for key,url in queries.items():
        rec={"url":url}
        try:
            params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"votable","QUERY":"SELECT TOP 3 * FROM rassmaster"} if key=="rassmaster" else {}
            r=requests.get(url,params=params,timeout=50)
            r.raise_for_status()
            if key=="rassmaster":
                import io
                t=Table.read(io.BytesIO(r.content),format="votable")
                rec["rows"]=[{k:str(row[k]) for k in t.colnames} for row in t]
            else:
                rec["text"]=r.text[:24000]
            rec["sha256"]=hashlib.sha256(r.content).hexdigest()
        except Exception as exc:
            rec["error"]=str(exc)
        (DOC/(key+".json")).write_text(json.dumps(rec,indent=2)+"\n")


def papers():
    terms=("XMM-Newton slew survey highly variable sources ROSAT Li Strotjohann",
           "eROSITA XMM slew transients XMMSL2 eRASS1 archival variability sources",
           "eROSITA ROSAT XMM long term extremely variable X-ray sources eRASS1 tidal disruption")
    for i,term in enumerate(terms):
        rec={"query":term}
        try:
            r=requests.get("https://api.crossref.org/works",params={"query.bibliographic":term,"rows":15,"select":"title,DOI,author,published"},timeout=45)
            r.raise_for_status()
            rec["works"]=[dict(title=x.get("title",[]),doi=x.get("DOI"),first_author=(x.get("author") or [{}])[0].get("family"),date=x.get("published")) for x in r.json()["message"]["items"]]
            rec["sha256"]=hashlib.sha256(r.content).hexdigest()
        except Exception as exc:
            rec["error"]=str(exc)
        (DOC/f"literature_{i}.json").write_text(json.dumps(rec,indent=2)+"\n")


if __name__=="__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        tasks = [pool.submit(inspect_source,item) for item in URLS.items()]
        tasks += [pool.submit(metadata),pool.submit(papers)]
        for t in tasks:
            t.result()
