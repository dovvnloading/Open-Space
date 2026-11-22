"""Test public archive routes on a common control before interpreting empty cones."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
from pathlib import Path
import tarfile

import requests
from astropy.io import fits
from astropy.table import Table

OUT=Path("results/xmmsl2-erass1-rass/outputs/service_audit")
OUT.mkdir(exist_ok=True,parents=True)
mpe="https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/"
args=[]

def query(label,url,params=None,method="GET",data=None):
    args.append((label,url,params or {},method,data))

for host in ("simbad.cds.unistra.fr","simbad.u-strasbg.fr","simbad.cfa.harvard.edu","simbad.harvard.edu"):
    tap=f"https://{host}/simbad/sim-tap/sync"
    q={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":"SELECT TOP 1 main_id,ra,dec FROM basic WHERE main_id='Vega'"}
    query("simbad_get_"+host,tap,q)
    query("simbad_post_"+host,tap,method="POST",data=q)
    query("simbad_cone_"+host,"https://"+host+"/simbad/sim-coo",{"Coord":"279.23473 +38.78369","Radius":"20","Radius.unit":"arcsec","output.format":"VOTable"})
for host in ("vizier.cfa.harvard.edu","vizier.cds.unistra.fr","vizier.nao.ac.jp","vizier.iucaa.in","vizier.inasan.ru"):
    query("vsx_"+host,f"https://{host}/viz-bin/asu-tsv",{"-source":"B/vsx/vsx","-c":"19 25 27.9 +42 47 03.7","-c.rs":20,"-out":"Name,RAJ2000,DEJ2000,Type","-out.max":8})
query("nodl_schema","https://datalab.noirlab.edu/query/query",{"sql":"select table_schema,table_name from information_schema.tables where lower(table_schema) like '%wise%' or lower(table_schema) like '%gaia%' or lower(table_schema) like '%tns%' order by table_name","ofmt":"csv","async":"False"})
for schema,fields in (("gaia_dr3.gaia_source","source_id,ra,dec,parallax,pmra,pmdec"),("gaia_edr3.gaia_source","source_id,ra,dec"),("allwise.source","designation,ra,dec"),("allwise.allwise_p3as_psd","designation,ra,dec"),("unwise_dr1.object","ra,dec")):
    query("datalab_"+schema,"https://datalab.noirlab.edu/query/query",{"sql":f"select {fields} from {schema} limit 1","ofmt":"csv","async":"False"})
for tag,url,table,cols in (("IRSA","https://irsa.ipac.caltech.edu/TAP/sync","allwise_p3as_psd","designation,ra,dec"),("Gaia","https://gea.esac.esa.int/tap-server/tap/sync","gaiadr3.gaia_source","source_id,ra,dec,parallax"),("GAIA_AIP","https://gaia.aip.de/tap/sync","gaiadr3.gaia_source","source_id,ra,dec"),("GAIA_ARI","https://gaia.ari.uni-heidelberg.de/tap/sync","gaiadr3.gaia_source","source_id,ra,dec")):
    query(tag,url,{"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":f"SELECT TOP 1 {cols} FROM {table}"})
query("tns_official_object","https://www.wis-tns.org/object/2018cow")
query("vsx_aavso","https://www.aavso.org/vsx/index.php",{"view":"api.list","ra":291.36625,"dec":42.78436,"radius":0.01,"format":"json"})
for slug in ("HamStar_eRASS1_Main_Likely_Identifications_v1.1","HamStar_eRASS1_Main_Possible_Identifications_v1.1"):
    query("model_"+slug,mpe+"FreundS_DR1/"+slug+".html")
    query("FITS_"+slug,mpe+"FreundS_DR1/"+slug+".fits.tgz")


def one(args):
    label,url,params,method,data=args
    rec={"label":label,"url":url,"params":params,"method":method,"data":data}
    try:
        headers={"User-Agent":"Open-Space public survey research"}
        if "datalab.noirlab.edu" in url:
            headers["X-DL-AuthToken"]="anonymous.0.0.anon_access"
        r=requests.request(method,url,params=params,data=data,headers=headers,timeout=(8,35))
        r.raise_for_status()
        rec.update(status=r.status_code,sha256=hashlib.sha256(r.content).hexdigest(),final_url=r.url)
        if label.startswith("FITS_"):
            with tarfile.open(fileobj=io.BytesIO(r.content)) as tar:
                f=next(m for m in tar if m.isfile() and m.name.lower().endswith((".fits",".fit")))
                with fits.open(io.BytesIO(tar.extractfile(f).read()),memmap=False) as tab:
                    h=tab[1]
                    rec["fields"]=h.columns.names
                    rec["example"]={k:str(h.data[k][0])[:80] for k in h.columns.names}
        else:
            rec["text"]=r.text[:8500]
            if b"<VOTABLE" in r.content[:1500] and b"QUERY_STATUS" in r.content and len(r.content)<1e6:
                try:
                    t=Table.read(io.BytesIO(r.content),format="votable")
                    rec["rows"]=len(t)
                    rec["fields"]=t.colnames
                    rec["example"]=str(t[0])[:1300] if len(t) else ""
                except Exception as e:
                    rec["vo_error"]=str(e)
    except Exception as e:
        rec["error"]=str(e)
    (OUT/(label+".json")).write_text(json.dumps(rec,indent=2)+"\n")
    print(label,rec.get("status"),rec.get("error"),flush=True)

with ThreadPoolExecutor(max_workers=7) as pool:
    list(pool.map(one,args))
