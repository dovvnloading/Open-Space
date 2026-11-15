"""Print archive queries for the frozen astrometry and fine ephemerides."""
import csv
from datetime import datetime,timedelta
from pathlib import Path
from urllib.parse import urlencode

root=Path(__file__).resolve().parent
sites={r['site']:r for r in csv.DictReader((root/'sites.csv').open())}
grids={('2027-02-06','16:47'):['IAR','Warkworth','McMurdo','SouthPole','GBT'],
       ('2054-03-02','04:00'):['MeerKAT','GMRT','Parkes','FAST']}
for (day,clock),keys in grids.items():
    start=datetime.fromisoformat(day+'T'+clock)
    end=start+timedelta(minutes=32)
    for key in keys:
        site=sites[key]
        lon,lat,h=[float(site[k]) for k in ('east_longitude_deg','latitude_deg','height_km')]
        jp=dict(COMMAND='299',CENTER='coord@399',MAKE_EPHEM='YES',EPHEM_TYPE='OBSERVER',
                START_TIME=str(start),STOP_TIME=str(end),STEP_SIZE='1 m',QUANTITIES='1,4,13,20,23',
                CSV_FORMAT='YES',ANG_FORMAT='DEG',EXTRA_PREC='YES',TIME_DIGITS='SECONDS',
                COORD_TYPE='GEODETIC',SITE_COORD=f'{lon},{lat},{h}')
        print(day,key,'JPL','https://ssd.jpl.nasa.gov/api/horizons.api?'+urlencode({'format':'json',**{k:f"'{v}'" for k,v in jp.items()}}))
        mi={'-name':'p:Venus','-ep':start.isoformat(),'-nbd':33,'-step':'1m','-tscale':'UTC',
            '-theory':'INPOP','-teph':1,'-tcoor':1,'-rplane':1,'-mime':'json','-observer':f'{lat},{lon},{h*1000}'}
        print(day,key,'INPOP','https://ssp.imcce.fr/webservices/miriade/api/ephemcc.php?'+urlencode(mi))
q="SELECT TOP 5 * FROM AS110.racs_low2_v02 WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',296.02183,-17.83641,0.01))"
print('RACS-low2','https://casda.csiro.au/casda_vo_tools/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q)))
