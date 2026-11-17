"""JPL historical known-object checks, with the actual WISE spacecraft state if available."""
import io,csv,json
from urllib.parse import urlencode
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from catalogs import HERE,query

def main():
    r=json.loads((HERE/'positions.json').read_text())[0]['visits']
    times=Time([x['mjd'] for x in r],format='mjd',scale='utc')
    url='https://ssd.jpl.nasa.gov/api/horizons.api?'+urlencode(dict(
        format='json',COMMAND="'-163'",CENTER="'500@399'",MAKE_EPHEM="'YES'",EPHEM_TYPE="'VECTORS'",
        TLIST="'"+','.join(str(v) for v in times.tdb.jd)+"'",OUT_UNITS="'KM-S'",REF_PLANE="'FRAME'",
        VEC_TABLE="'2'",CSV_FORMAT="'YES'",VEC_CORR="'NONE'"))
    query(('wise_horizons',url))
    result=json.loads((HERE/'wise_horizons.response').read_text())['result']
    states=[{} for v in r]
    if '$$SOE' in result:
        line=result.split('$$SOE')[1].split('$$EOE')[0].strip()
        states=[{'xobs':','.join(str(float(v)) for v in row[2:8])} for row in csv.reader(io.StringIO(line))]
    else:states=[{'mpc-code':'500'} for v in r]
    jobs=[]
    for v,row in enumerate(r):
        c=SkyCoord(row['ra']*u.deg,row['dec']*u.deg)
        url='https://ssd-api.jpl.nasa.gov/sb_ident.api?'+urlencode({
            **states[v],'obs-time':times[v].isot[:19],
            'fov-ra-center':c.ra.to_string(unit=u.hourangle,sep='-',precision=3,pad=True),
            'fov-dec-center':c.dec.to_string(unit=u.deg,sep='-',precision=2,pad=True),
            'fov-ra-hwidth':.09,'fov-dec-hwidth':.065,'mag-required':'false','two-pass':'true',
            'suppress-first-pass':'true','req-elem':'false'})
        jobs.append(('jpl'+str(v),url))
    for job in jobs:
        query(job)

if __name__=='__main__':main()
