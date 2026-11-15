"""Archive original PR #24 light curves and try independent, *unselected* holdouts.

The output is uploaded as a GitHub Actions artifact; failed external surveys
are recorded rather than silently interpreted as nondetections.
"""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile

RA, DEC = 277.941193, -14.315539
OUT = Path(os.environ.get("V612_HARVEST_OUT", "/tmp/v612-harvest"))
OUT.mkdir(parents=True,exist_ok=True)
status = {}
for sector in (80,92):
    query = urlencode(dict(ra=RA,dec=DEC,y=5,x=5,units="px",sector=sector))
    url = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + query
    try:
        with urlopen(url,timeout=175) as response:
            raw = response.read()
        with ZipFile(io.BytesIO(raw)) as z:
            name = next(n for n in z.namelist() if n.endswith(".fits"))
            content = z.read(name)
        (OUT / f"tess{sector}.fits").write_bytes(content)
        status[f"TESS {sector}"]={"status":"downloaded","sha256":hashlib.sha256(content).hexdigest(),
                                    "bytes":len(content),"query_url":url}
    except Exception as e:
        status[f"TESS {sector}"]={"status":"failed","error":repr(e),"query_url":url}

try:
    from pyasassn.client import SkyPatrolClient
    client = SkyPatrolClient(verbose=False)
    cat = client.cone_search(RA, DEC, radius=9, units="arcsec",catalog="master_list")
    cat.to_csv(OUT/"asas-catalog.csv",index=False)
    status["ASAS-SN cone"]={"status":"catalog_retrieved","n":len(cat)}
    cat = cat.copy()
    cat["offset_arcsec"] = ((cat["ra_deg"]-RA)*3490.0)**2+((cat["dec_deg"]-DEC)*3600.0)**2
    cat=cat.sort_values("offset_arcsec").head(4)
    for _,row in cat.iterrows():
        ident = int(row["asas_sn_id"])
        try:
            lc = client.query_list([ident],catalog="master_list",download=True)
            curve = lc.data
            if curve is None or len(curve)==0:
                raise ValueError("No photometric rows")
            name=f"asas-{ident}.csv"
            curve.to_csv(OUT/name,index=False)
            status[name] = {"status":"downloaded","n":len(curve),
                            "ra":float(row["ra_deg"]),"dec":float(row["dec_deg"])}
        except Exception as e:
            status[f"asas-{ident}"] = {"status":"failed","error":repr(e)}
except Exception as e:
    status["ASAS-SN cone"]={"status":"failed","error":repr(e)}

(OUT/"acquisition.json").write_text(json.dumps(status,indent=2)+"\n")
print("ACQUISITION_SUMMARY="+json.dumps(status,sort_keys=True),flush=True)
