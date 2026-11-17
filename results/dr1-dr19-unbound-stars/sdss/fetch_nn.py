"""Fetch parallax-free astroNN predictions for APOGEE screen sources; weighted_dist uses Gaia."""
from concurrent.futures import ThreadPoolExecutor
import csv
from pathlib import Path
from fetch_survey import HERE,query
def get(ids):
 q=("SELECT sdss_id,gaia_dr3_source_id,release,telescope,spectrum_pk,nn_dist,nn_dist_err,"
    "nn_dist_model_err,weighted_dist,weighted_dist_err,result_flags,flag_warn,flag_bad "
    "FROM astro_nn_dist_apogee_star WHERE release='sdss5' AND sdss_id IN ("+ids+") ORDER BY sdss_id,spectrum_pk")
 return query(q).splitlines()[1:]
if __name__=='__main__':
 r=list(csv.DictReader((HERE/'prefilter_tail.csv').open()))
 # Include independently RV-supported targets near the follow-up speed threshold.
 ids=sorted({x['sdss_id'] for x in r if 'apogee' in x['survey']}|{'78976052','88803696','86480083','86409747'},key=int)
 chunks=[','.join(ids[i:i+60]) for i in range(0,len(ids),60)]
 with ThreadPoolExecutor(5) as p: chunks=list(p.map(get,chunks))
 buf=[chunks[0][0]]+[v for c in chunks for v in c[1:]]
 (HERE/'astro_nn_apogee.csv').write_bytes(b'\n'.join(buf)+b'\n')
 print(len(buf)-1)
