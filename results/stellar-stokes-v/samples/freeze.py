"""Build four epoch-aware samples and preserve per-planet periods from public releases."""
import argparse
import csv
import gzip
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd

HERE = Path(__file__).resolve().parent
GAIA = 'https://gea.esac.esa.int/tap-server/tap/sync'
NASA = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync'
CNS5 = 'https://dc.g-vo.org/tap/sync'
CENSUS = 'https://exoplanetarchive.ipac.caltech.edu/data/Contributed/20pcCensus/Kirkpatrick2024_20pcCensus_Table4_Parameters.csv'
NEAR = 'select source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,parallax_over_error,ruwe from gaiadr3.gaia_source where parallax>=20 and parallax_over_error>=10 and pmra is not null and pmdec is not null'
PLANETS = 'select pl_name,hostname,ra,dec,gaia_dr3_id,sy_dist,pl_orbper,pl_orbpererr1,pl_tranmid,pl_orbtper,st_rotp,sy_pmra,sy_pmdec from pscomppars order by hostname,pl_name'
CNS = 'select cns5_id,gaia_edr3_id,ra,dec,epoch,pmra,pmdec,parallax from cns5.main order by cns5_id'
FIELDS = ['hostname','gaia_dr3_id','ra_deg','dec_deg','ref_epoch','pmra_masyr','pmdec_masyr','dist_pc','astrometry']
NAMES = dict(ra='ra_deg',dec='dec_deg',pmra='pmra_masyr',pmdec='pmdec_masyr',source_id='gaia_dr3_id')
def tap(base,query,maximum):
    return base+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=maximum,QUERY=query))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cache',type=Path)
    ap.add_argument('out',type=Path,help='plain CSVs suitable for the two catalogue matchers')
    args=ap.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    args.out.mkdir(parents=True,exist_ok=True)
    receipts={}
    def fetch(item):
        name,url=item
        path=args.cache/(name+'.csv')
        if not path.exists():
            with urlopen(url,timeout=240) as response:
                raw=response.read()
            if b',' not in raw.splitlines()[0] or raw.lstrip().startswith(b'<'):
                raise RuntimeError('Invalid catalogue response '+name)
            path.write_bytes(raw)
        raw=path.read_bytes()
        receipts[name]=dict(url=url,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
        return pd.read_csv(path,dtype=str,keep_default_na=False)
    planets,near,cns,census=list(map(fetch,[
        ('hosts_planets',tap(NASA,PLANETS,150000)),('gaia_nearby',tap(GAIA,NEAR,150000)),
        ('cns5',tap(CNS5,CNS,15000)),('census',CENSUS)]))
    planets['gaia_dr3_id']=planets.gaia_dr3_id.str.replace('Gaia DR3 ','',regex=False)
    hosts=planets.drop_duplicates('hostname').copy()
    ids=sorted(set(planets.gaia_dr3_id)-{''})
    jobs=[]
    for i in range(0,len(ids),170):
        query='select source_id,ra,dec,ref_epoch,pmra,pmdec,parallax from gaiadr3.gaia_source where source_id in ('+','.join(ids[i:i+170])+') order by source_id'
        jobs.append((f'gaia_host_{i:04d}',tap(GAIA,query,1000)))
    with ThreadPoolExecutor(4) as pool:
        gaia=pd.concat(list(pool.map(fetch,jobs))).set_index('source_id',verify_integrity=True)
    hosts=hosts.rename(columns={'sy_dist':'dist_pc','ra':'archive_ra','dec':'archive_dec'})
    hosts=hosts.join(gaia,on='gaia_dr3_id')
    usable=hosts.pmra.fillna('').ne('') & hosts.pmdec.fillna('').ne('')
    hosts['astrometry']='Gaia DR3'
    hosts.loc[~usable,'astrometry']='NASA archive position; epoch or motion unverified'
    for key in ['ra','dec']:
        hosts.loc[~usable,key]=hosts.loc[~usable,'archive_'+key]
    hosts.loc[~usable,['ref_epoch','pmra','pmdec']]=''
    near=near.sort_values('source_id')
    near['hostname']='Gaia DR3 '+near.source_id
    near['astrometry']='Gaia DR3'
    near['dist_pc']=1000/pd.to_numeric(near.parallax)
    cns=cns[cns.ra.ne('') & cns.dec.ne('')].copy()
    cns['hostname']='CNS5 '+cns.cns5_id
    cns['astrometry']='CNS5 (native epoch)'
    cns['dist_pc']=1000/pd.to_numeric(cns.parallax)
    cns=cns.rename(columns={'epoch':'ref_epoch','gaia_edr3_id':'source_id'})
    # L0 through Y. Exclude white dwarfs (indices above 100) and system rows
    # without their own measured coordinates. A spectral type is not a mass.
    incl=pd.Series(False,index=census.index)
    for col in ['SpTOpt_indx','SpTNIR_indx']:
        incl |= pd.to_numeric(census[col],errors='coerce').between(50,79.99)
    d=census[incl].rename(columns={'DefaultName':'hostname','RA (deg)':'ra','Dec (deg)':'dec','Epoch':'ref_epoch','PMRA (mas/yr)':'pmra','PMDec (mas/yr)':'pmdec'}).copy()
    miss=d[['ra','dec','ref_epoch','pmra','pmdec']].eq('').any(axis=1)
    d=d[~miss].copy()
    d['dist_pc']=1000/pd.to_numeric(d['Parallax (mas)'],errors='coerce')
    d['source_id']=d.Gaia.str.extract(r'Gaia (?:E?DR3) (\d+)',expand=False)
    d['astrometry']='Kirkpatrick 2024 Table 4 (native epoch)'
    products={}
    for label,frame in [('hosts',hosts),('gaia_nearby',near),('cns5',cns),('dwarfs',d)]:
        frame=frame.rename(columns=NAMES)[FIELDS].sort_values('hostname')
        data=frame.to_csv(index=False,lineterminator='\n').encode()
        (args.out/(label+'.csv')).write_bytes(data)
        products[label]=dict(rows=len(frame),sha256=hashlib.sha256(data).hexdigest())
        if label in ['hosts','dwarfs']:
            (HERE/(label+'.csv.gz')).write_bytes(gzip.compress(data,mtime=0))
    pcols=['pl_name','hostname','pl_orbper','pl_orbpererr1','pl_tranmid','pl_orbtper','st_rotp']
    data=planets[pcols].to_csv(index=False,lineterminator='\n').encode()
    (HERE/'periods.csv.gz').write_bytes(gzip.compress(data,mtime=0))
    doc=dict(inputs=receipts,outputs=products,planets=len(planets),hosts_with_motion=int(usable.sum()),
             spectral_LTY_before_astrometry_filter=int(incl.sum()),dwarfs_missing_astrometry=int(miss.sum()))
    (HERE/'receipt.json').write_text(json.dumps(doc,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':
    main()
