"""JPL two-pass check at the middle of the four short Euclid pointings."""
import csv,io,json
from urllib.parse import urlencode
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from common import HERE,RA0,DE0,T
from metadata import retrieve

def main():
    t=Time(T.mjd.mean(),format='mjd',scale='utc')
    url='https://ssd.jpl.nasa.gov/api/horizons.api?'+urlencode(dict(
        format='json',COMMAND="'-680'",CENTER="'500@399'",MAKE_EPHEM="'YES'",EPHEM_TYPE="'VECTORS'",TLIST=f"'{t.tdb.jd}'",OUT_UNITS="'KM-S'",REF_PLANE="'FRAME'",VEC_TABLE="'2'",CSV_FORMAT="'YES'",VEC_CORR="'NONE'"))
    result=json.loads(retrieve(('jpl_observer',url)).read_bytes())['result'].split('$$SOE')[1].split('$$EOE')[0].strip()
    state=next(csv.reader(io.StringIO(result)))[2:8]
    c=SkyCoord(RA0*u.deg,DE0*u.deg)
    q={'obs-time':t.isot[:19],'xobs':','.join(str(float(i)) for i in state),'fov-ra-center':c.ra.to_string(unit=u.hourangle,sep='-',precision=3,pad=True),'fov-dec-center':c.dec.to_string(unit=u.deg,sep='-',precision=2,pad=True),'fov-ra-hwidth':.025,'fov-dec-hwidth':.025,'two-pass':'true','mag-required':'false','suppress-first-pass':'true','req-elem':'false'}
    retrieve(('small_bodies','https://ssd-api.jpl.nasa.gov/sb_ident.api?'+urlencode(q)))
if __name__=='__main__':main()
