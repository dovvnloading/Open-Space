"""Check the official Q1 SPE galaxy redshift and fitted emission-line tables by exact MER ID.

python spe_check.py MODEL_CATALOG_DIR RESPONSE_CACHE_DIR
"""
import csv
import hashlib
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from infer import HERE, read, write

ENDPOINT='https://eas.esac.esa.int/tap-server/tap/sync'
GAL='q1.spectro_zcatalog_spe_galaxy_candidates'
LINES='q1.spectro_line_features_catalog_spe_line_features_cat'
ZCOL=['object_id','spe_rank','spe_z','spe_z_err','spe_z_prob','spe_z_rel','spe_subclass','spe_cont_snr']
LCOL=['object_id','spe_rank','spe_line_id','spe_line_name','spe_line_central_wl_gf','spe_line_flux_gf',
      'spe_line_snr_gf','spe_line_flux_di','spe_line_snr_di','spe_line_flag','spe_line_n_dith']


def run(data,cache):
    cache.mkdir(parents=True,exist_ok=True)
    cat = read(data/'q1_discovery_engine_lens_catalog.csv')
    primary = {r['object_id'] for r in cat if r['subset']=='discovery_engine' and r['grade'] in ('A','B')}
    models = {r['id_str'].split('_')[-1].replace('NEG','-') for r in read(data/'modeling_lens_mass.csv')}
    cone_ids = {r['object_id'] for r in read(HERE/'cone_matches.csv') if r['survey']=='euclid'}
    ids=sorted(primary|models|cone_ids)
    def query(job):
        tab,cols,condition,keys=job
        sql='SELECT '+','.join(cols)+' FROM '+tab+' WHERE spe_rank=0 AND object_id IN ('+','.join(keys)+')'+condition+' ORDER BY object_id'
        query_sha=hashlib.sha256(sql.encode()).hexdigest()
        f=cache/(query_sha+'.csv')
        if not f.exists():
            url=ENDPOINT+'?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql))
            for attempt in range(3):
                try:
                    b=urllib.request.urlopen(url,timeout=100).read()
                    if not b.startswith(b'object_id,'):
                        raise ValueError(b[:350])
                    f.write_bytes(b)
                    break
                except (OSError,ValueError):
                    if attempt==2:
                        raise
                    time.sleep(1)
        raw=f.read_bytes()
        return dict(table=tab,query=sql,query_sha256=query_sha,response_sha256=hashlib.sha256(raw).hexdigest(),
                    n_ids=len(keys),n_rows=len(read(f))),read(f)
    def batches(tab,cols,cond,keys):
        return [(tab,cols,cond,keys[i:i+45]) for i in range(0,len(keys),45)]
    with ThreadPoolExecutor(4) as pool:
        zd=list(pool.map(query,batches(GAL,ZCOL,'',ids)))
        z_rows=[r for _,rows in zd for r in rows]
        special = cone_ids|{r['object_id'] for r in z_rows if .9<float(r['spe_z'])<1.8 and float(r['spe_z_prob'])>=.99}
        line_ids=sorted(special & {r['object_id'] for r in z_rows})
        ld=list(pool.map(query,batches(LINES,LCOL,'',line_ids)))
    ls=[r for _,rows in ld for r in rows]
    def numeric(v):
        return float(v) if v else float('nan')
    successes=set()
    for l in ls:
        name=l['spe_line_name'].replace('-','').replace('_','').lower()
        wl=numeric(l['spe_line_central_wl_gf'])
        l['is_Halpha']=int(name in ('halpha','ha','hα'))
        l['Halpha_GF_pass']=int(l['is_Halpha'] and numeric(l['spe_line_flux_gf'])>=2e-16 and
                                numeric(l['spe_line_snr_gf'])>=5 and (12500<=wl<=18500))
        if l['Halpha_GF_pass']:
            successes.add(l['object_id'])
    for z in z_rows:
        z['in_Q1_primary']=int(z['object_id'] in primary)
        z['in_cones']=int(z['object_id'] in cone_ids)
        z['pass_FAQ']=int(.9<numeric(z['spe_z'])<1.8 and numeric(z['spe_z_prob'])>=.99 and z['object_id'] in successes)
    write(HERE/'spe_rank0.csv',z_rows)
    write(HERE/'spe_lines.csv',ls)
    (HERE/'spe_manifest.json').write_text(json.dumps(dict(n_primary=len(primary),n_model_ids=len(models),
      n_cone_ids=len(cone_ids),n_distinct_queries=len(ids),n_rank0=len(z_rows),n_FAQ_pass=sum(x['pass_FAQ'] for x in z_rows),
      queries=[m for m,_ in zd+ld]),indent=1)+'\n')


if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]))
