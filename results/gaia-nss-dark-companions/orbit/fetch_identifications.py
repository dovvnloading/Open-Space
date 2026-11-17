"""Resolve source identities and the HD 264291 bibliography through SIMBAD TAP."""
import hashlib
import json
from urllib.parse import urlencode
from urllib.request import urlopen
from pathlib import Path

HERE=Path(__file__).resolve().parent
TAP="https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?"
QUERIES={
 "candidate_identifications.csv":"""SELECT a.id,b.main_id,b.otype,b.sp_type,b.ra,b.dec FROM ident AS a JOIN basic AS b ON b.oid=a.oidref WHERE id IN ('Gaia DR3 3378588057203660160','Gaia DR3 1830677037678583680','Gaia DR3 4661290764764683776','Gaia DR3 275582835635329024','Gaia DR3 902085660651149440','Gaia DR3 2077800002648034944','Gaia DR3 2130638481326861824','Gaia DR3 4658457941738020736','Gaia DR3 5713565792393224960')""",
 "hd264291_bibliography.csv":"""SELECT r.bibcode,r.title FROM ident AS a JOIN has_ref AS h ON a.oidref=h.oidref JOIN ref AS r ON r.oidbib=h.oidbibref WHERE a.id='Gaia DR3 3378588057203660160'"""
}
if __name__=="__main__":
    hashes={}
    for file,query in QUERIES.items():
        url=TAP+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query))
        with urlopen(url,timeout=40) as response:
            raw=response.read()
        if b"HD 264291" not in raw and b"bibcode,title" not in raw:
            raise ValueError("Unexpected SIMBAD response")
        (HERE/file).write_bytes(raw)
        hashes[file]=dict(url=url,sha256=hashlib.sha256(raw).hexdigest())
    (HERE/"identification_inputs.json").write_text(json.dumps(hashes,indent=2)+"\n")
