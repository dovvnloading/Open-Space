#!/usr/bin/env python3
"""Bounded access checks for external Gaia, SIMBAD and alternate PPS archive.
No science claim can be inferred from a service probe.
"""
import json
from pathlib import Path
import requests
import time
import datetime as dt

BASE=Path(__file__).resolve().parent
OUT=BASE/"outputs"
S=requests.Session()
S.headers.update({"User-Agent":"Open-Space-QPE-survey-service-audit/1.0"})

def trial(name,url,params=None,method="GET"):
    t0=time.monotonic()
    out={"name":name,"url":url}
    try:
        with S.request(method,url,params=params,stream=True,timeout=(9,19),allow_redirects=True) as r:
            data=next(r.iter_content(2048),b"")
            out.update(http_status=r.status_code,content_type=r.headers.get("Content-Type"),
                    content_length=r.headers.get("Content-Length"),
                    head=data[:1000].decode("utf-8","replace"),
                    bytes_first_block=len(data),final_url=r.url)
    except Exception as e:out["error"]=repr(e)[:250]
    out["elapsed_s"]=round(time.monotonic()-t0,2)
    return out

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg.get("stage")=="archive_services"
    gaia=("SELECT TOP 5 source_id,ra,dec,parallax,parallax_over_error,pmra,pmdec "
          "FROM gaiadr3.gaia_source WHERE 1=CONTAINS("
          "POINT('ICRS',ra,dec),CIRCLE('ICRS',195.50054,27.78267,0.00056))")
    simbad=("SELECT TOP 5 main_id,ra,dec,otype FROM basic WHERE 1=CONTAINS("
            "POINT('ICRS',ra,dec),CIRCLE('ICRS',195.50054,27.78267,0.00056))")
    tap={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":gaia}
    results=[
        trial("esa_gaia", "https://gea.esac.esa.int/tap-server/tap/sync",tap),
        trial("ari_gaia", "https://gaia.ari.uni-heidelberg.de/tap/sync",tap),
        trial("vizier_gaia", "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
              {**tap,"QUERY":"SELECT TOP 5 Source,RA_ICRS,DE_ICRS,Plx,pmRA,pmDE FROM \"I/355/gaiadr3\" WHERE 1=CONTAINS(POINT('ICRS',RA_ICRS,DE_ICRS),CIRCLE('ICRS',195.50054,27.78267,0.00056))"}),
        trial("simbad", "https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
              {**tap,"QUERY":simbad}),
        trial("heasarc_xmm_directory","https://heasarc.gsfc.nasa.gov/FTP/xmm/data/rev0/"),
        trial("alternate_xsa", "https://xsa.esac.esa.int/xsa-sl/servlet/data-action-aio",
              {"obsno":"0872392901","level":"PPS","instname":"PN",
               "name":"SRCTSR","extension":"FTZ"})
    ]
    report={"utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "status":"SERVICE_PROBES_ONLY","results":results,
            "archive_screen_permitted":False}
    OUT.joinpath("archive_service_probes.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps([{k:r.get(k) for k in ("name","http_status","elapsed_s","error")}
         for r in results],indent=2),flush=True)

if __name__=="__main__":main()
