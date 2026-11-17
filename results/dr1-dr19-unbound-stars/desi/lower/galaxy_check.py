from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from urllib.parse import urlencode
import csv
from pathlib import Path
from math import cos,radians,atan2,sin,sqrt,asin,degrees
import pandas as pd
ROOT=Path(__file__).resolve().parent
leads=pd.read_csv(ROOT/'lower_all.csv',dtype={'source_id':str});lead=leads[leads.source_id.isin(['1070549357629086976','1525802834304705024','1545579445011222272','1717410430709663616','4444095982221534720'])]
fields='main_id,ra,dec,otype,rvz_radvel,galdim_majaxis,galdim_minaxis,galdim_angle'
u='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
def get(r):
 q=f"SELECT {fields} FROM basic WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r.ra},{r.dec},0.25)) AND galdim_majaxis>1"
 url=u+'?'+urlencode({'request':'doQuery','lang':'adql','format':'csv','query':q})
 try:b=urlopen(url,timeout=60).read();print(r.source_id,b[:270],flush=True);(ROOT/('galaxies_'+r.source_id+'.csv')).write_bytes(b)
 except Exception as e:print(str(e))
if __name__=='__main__':
    with ThreadPoolExecutor(4) as pool: list(pool.map(get,lead.itertuples()))
    q="select main_id,ra,dec,rvz_radvel,galdim_majaxis,galdim_minaxis from basic join ident on oid=oidref where id in ('M 82','M 81','M 63','M 106')"
    url=u+'?'+urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q})
    with urlopen(url,timeout=60) as f: (ROOT/'bright_galaxies.csv').write_bytes(f.read())

