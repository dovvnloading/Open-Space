"""Query public exposure indices around all four modeled event times."""
import csv
import io
import json
import math
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from urllib.parse import urlencode
from urllib.request import Request, urlopen

FIELDS = {
  '2024': ('2024-02-07T17:49:04',356.4842526,-2.8498067),
  '2026': ('2026-01-31T00:34:52',0.2954826,-1.3025397),
}
TAP = {'ESO':'https://archive.eso.org/tap_obs/sync',
       'CADC':'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync',
       'MAST':'https://mast.stsci.edu/vo-tap/api/v0.1/caom/sync'}

def search(task):
    year, archive = task
    date,ra,dec=FIELDS[year]
    t=datetime.fromisoformat(date).replace(tzinfo=timezone.utc).timestamp()
    mjd=t/86400+40587
    # +/- 60 min brackets the full four-model range on both nights.
    if archive in TAP:
        query=(f"SELECT TOP 200 obs_collection,obs_id,t_min,t_max FROM ivoa.ObsCore "
               f"WHERE 1=INTERSECTS(s_region,CIRCLE('ICRS',{ra},{dec},0.003)) "
               f"AND t_min<{mjd+1/24:.8f} AND t_max>{mjd-1/24:.8f}")
        args=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)
        url=TAP[archive]+'?'+urlencode(args)
        req=url
    else:
        url='https://astroarchive.noirlab.edu/api/adv_search/find/?rectype=file'
        labels=[datetime.fromtimestamp(t+h,tz=timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
                for h in (-3600,3600)]
        args=dict(outfields=['archive_filename','dateobs_center','ra_center','dec_center','instrument'],
                  search=[['dateobs_center',*labels],['dec_center',dec-4,dec+4]])
        req=Request(url,data=json.dumps(args).encode(),headers={'Content-Type':'application/json'})
    result=dict(year=year,archive=archive,url=url,request=args)
    try:
        raw=urlopen(req,timeout=80).read().decode()
        if archive in TAP:
            header=next(csv.reader(io.StringIO(raw)))
            if 'obs_id' not in header:raise ValueError(raw[:250])
            rows=list(csv.DictReader(io.StringIO(raw)))
            result['matches']=rows
        else:
            rows=json.loads(raw)[1:]
            result['declination_band_rows']=len(rows)
            # Generous centre-distance fallback for the ground-based fields.
            def angle(r):
                a=math.radians(float(r['dec_center']));b=math.radians(dec)
                x=math.radians(float(r['ra_center'])-ra)
                return math.degrees(math.acos(min(1,max(-1,math.sin(a)*math.sin(b)+math.cos(a)*math.cos(b)*math.cos(x)))))
            result['matches']=[r for r in rows if angle(r)<4]
            if rows:result['nearest_boresight_deg']=round(min(map(angle,rows)),3)
    except Exception as e:
        result['error']=str(e)
    return result

if __name__=='__main__':
    tasks=[(y,a) for y in FIELDS for a in [*TAP,'NOIRLab']]
    with ThreadPoolExecutor(max_workers=4) as ex:
        print(json.dumps(list(ex.map(search,tasks)),indent=2))
