"""Save unchanged Hydrocron CSV extracts, paired by reach and collection; HTTP metadata in manifest."""
from pathlib import Path
import concurrent.futures, csv, datetime, hashlib, io, json, sys, time
from urllib.parse import urlencode
from urllib.request import urlopen, Request

out=Path(__file__).resolve().parent
out.mkdir(exist_ok=True)
reaches="14277700101 14277700111 14277700091 14279000471 14291000011 14281000011 14293000101 14293000111 14293000091".split()
fields='reach_id,time_str,range_start_time,p_lon,p_lat,river_name,reach_q,reach_q_b,wse,wse_u,width,width_u,slope,slope_u,slope2,slope2_u,n_good_nod,obs_frac_n,partial_f,dark_frac,node_dist,loc_offset,xtrk_dist,cycle_id,pass_id,crid,sword_version,granuleUR'
node_fields='node_id,reach_id,time_str,range_start_time,lat,lon,node_q,node_q_b,wse,wse_u,width,width_u,n_good_pix,dark_frac,node_dist,xtrk_dist,cycle_id,pass_id,crid,sword_version,granuleUR'
url='https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries'
def one(job):
    feature,ID,version=job
    name=f'{feature.lower()}_{ID}_{version}.csv'
    params={'feature':feature,'feature_id':ID,'collection_name':'SWOT_L2_HR_RiverSP_'+version,'start_time':'2023-07-21T00:00:00Z','end_time':'2026-11-16T00:00:00Z','output':'csv','fields':fields if feature=='Reach' else node_fields}
    if feature=='Reach' and version=='D':
        cross=json.loads((out/'sword_crosswalk.json').read_text())['selected']
        params['feature_id']=next(r['v17_reach_id'] for r in cross if r['v16_reach_id']==ID)
    if feature=='Node' and version=='D':
        cross=json.loads((out/'sword_node_crosswalk.json').read_text())['selected']
        params['feature_id']=next(r['v17_node_id'] for r in cross if r['v16_node_id']==ID)
    del params['output']
    target=url+'?'+urlencode(params)
    for attempt in range(3):
        try:
            with urlopen(Request(target,headers={'Accept':'text/csv'}),timeout=70) as response:
                raw=response.read(); requested_date=response.headers.get('date'); rqid=response.headers.get('x-amzn-requestid')
            if raw.startswith(b'{'): raise ValueError(raw[:250])
            transport_sha256=hashlib.sha256(raw).hexdigest()
            if raw.startswith(b'"'): raw=json.loads(raw).encode('utf-8')
            rows=list(csv.DictReader(io.StringIO(raw.decode())))
            if not rows or rows[0].get('reach_id') is None: raise ValueError(raw[:250])
            (out/name).write_bytes(raw)
            print(name,len(rows),flush=True)
            return dict(file=name,http_response_sha256=transport_sha256,feature_id_queried=params['feature_id'],source_url=target,accept='text/csv',rows=len(rows),sha256=hashlib.sha256(raw).hexdigest(),http_date=requested_date,request_id=rqid,units='provider unit columns in CSV')
        except Exception as exc:
            print(name,repr(exc),flush=True)
            if attempt==2: raise
            time.sleep(1+attempt)
if __name__=='__main__':
    node_ids=sys.argv[1:]
    work=[('Node',ID,v) for ID in node_ids for v in ['D','2.0']] if node_ids else [('Reach',ID,v) for ID in reaches for v in ['D','2.0']]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool: manifest=list(pool.map(one,work))
    (out/('node_manifest.json' if node_ids else 'geometry_manifest.json')).write_text(json.dumps(manifest,indent=2)+'\n')
