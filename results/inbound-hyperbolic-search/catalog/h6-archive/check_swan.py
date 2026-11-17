"""Compare archived SWAN positions with Horizons and two published MPC orbit fits."""
import argparse
import csv
from datetime import date
import hashlib
import io
import json
import math
from pathlib import Path
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen

parser = argparse.ArgumentParser()
parser.add_argument("--cache", type=Path, help="directory of frozen or previously downloaded responses")
parser.add_argument("--output", type=Path, default=Path(__file__).with_name("swan-receipts.json"))
args = parser.parse_args()
receipts = []
archive = "https://soho.nascom.nasa.gov/data/summary/swan/comets/"
api = "https://ssd.jpl.nasa.gov/api/horizons.api?"

def read(name, url):
    path = args.cache / name if args.cache else None
    if path and path.exists():
        raw = path.read_bytes()
    else:
        request = Request(url, headers={"User-Agent": "Open-Space-SWAN-check/1.0"})
        with urlopen(request, timeout=60) as response:
            raw = response.read()
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
    receipts.append({"input": name, "url": url, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()})
    return raw

def mpc(code):
    return "https://www.minorplanetcenter.net/iau/mpec/K04/K04" + code + ".html"

def table(raw):
    body = json.loads(raw)["result"].split("$$SOE")[1].split("$$EOE")[0].strip()
    return list(csv.reader(io.StringIO(body), skipinitialspace=True))

def sphsep(a,b):
    cross = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    return math.degrees(math.atan2(math.hypot(*cross), sum(x*y for x,y in zip(a,b))))

def unit(ra,dec):
    a,d = map(math.radians,(ra,dec))
    return [math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]

def geodesic(row, ra_deg, dec_deg):
    h,m,s = map(float,row[3].split())
    ra = 15*(h+m/60+s/3600)
    deg,m,s = map(float,row[4].split())
    dec = math.copysign(abs(deg)+m/60+s/3600, deg)
    return sphsep(unit(ra,dec),unit(ra_deg,dec_deg))

def direction(t, observer, elements):
    q,e,tp,arg,node,inc = elements
    w,o,i = map(math.radians,(arg,node,inc))
    # Barker and hyperbolic Kepler equations in a common ecliptic J2000 frame.
    mu = 0.0002959122082855911
    dt = 0.0
    for _ in range(6):
        n = t - dt - tp
        if e == 1:
            d = 2*math.sinh(math.asinh(1.5*n*math.sqrt(mu/(2*q**3)))/3)
            x,y = q*(1-d*d),2*q*d
        else:
            a = q/(e-1)
            mean = n*math.sqrt(mu/a**3)
            h = mean/(e-1)
            for _ in range(20):
                step = (e*math.sinh(h)-h-mean)/(e*math.cosh(h)-1)
                h -= step
                if abs(step)<1e-13: break
            x,y = a*(e-math.cosh(h)),a*math.sqrt(e*e-1)*math.sinh(h)
        u,v = math.cos(w)*x-math.sin(w)*y, math.sin(w)*x+math.cos(w)*y
        heliocentric = (math.cos(o)*u-math.sin(o)*math.cos(i)*v,
                        math.sin(o)*u+math.cos(o)*math.cos(i)*v, math.sin(i)*v)
        los = [a-b for a,b in zip(heliocentric,observer)]
        dt = math.hypot(*los)/173.144632674240
    return los

iauc = read("h6-src0.cache","https://tamkin3.eps.harvard.edu/iauc/08300/08346.html")
june = read("h6-mpecl.cache",mpc("L42"))
july = read("h6-mpec.cache",mpc("N33"))
assert b"May 13" in iauc and b"0.77605" in june and b"-0.000124" in july
index = read("h6-swan.cache",archive).decode()
filenames = sorted(set(re.findall(r'swan_ARecl_sm_2004(?:04(?:2[89]|30)|05(?:0[1-9]|1[0-7]))_[0-9]+\.gif',index)))
for name in filenames:
    image = read("h6sky/"+name,archive+name)
    assert image[:3]==b"GIF" and int.from_bytes(image[6:8],"little")==640
obs_params = dict(format="json", COMMAND="'DES=2004 H6;'", OBJ_DATA="'YES'", MAKE_EPHEM="'YES'",
                  EPHEM_TYPE="'OBSERVER'", CENTER="'500@-21'", START_TIME="'2004-04-28'",
                  STOP_TIME="'2004-05-15'",STEP_SIZE="'1d'",QUANTITIES="'1,3,36'",CSV_FORMAT="'YES'")
vec_params = dict(format="json", COMMAND="'-21'", OBJ_DATA="'NO'", MAKE_EPHEM="'YES'",
                  EPHEM_TYPE="'VECTORS'", CENTER="'500@10'", START_TIME="'2004-04-28'",
                  STOP_TIME="'2004-05-15'",STEP_SIZE="'1d'",VEC_TABLE="'1'",REF_PLANE="'ECLIPTIC'",
                  OUT_UNITS="'AU-D'",CSV_FORMAT="'YES'")
ephemeris = table(read("h6-horizons.cache",api+urlencode(obs_params)))
vectors = table(read("h6-obsvec.cache",api+urlencode(vec_params)))
# The circular reports only an integer day, hour and minute of RA, and decimal degrees of Dec.
published = {"Apr-29":(4,9,-9.4),"May-02":(4,5,-10.3),"May-05":(3,59,-10.6),"May-08":(3,56,-11.4)}
may12 = date(2004,5,12).toordinal()+1721424.5
models = [(0.77605,1.0,may12+.680,269.030,317.239,107.653),
          (0.7762059,1.0004976,may12+.68449,269.04100,317.24449,107.66093)]
results=[]
for obs,vec in zip(ephemeris,vectors):
    d = obs[0].strip()[5:11]
    if d not in published: continue
    h,m,de = published[d]
    jd = float(vec[0]); observer = list(map(float,vec[2:5]))
    sep = sphsep(*(direction(jd,observer,model) for model in models))*3600
    results.append({"day": d, "IAUC8346_RA_hours_minutes":[h,m],"IAUC8346_Dec_deg":de,
                    "SWAN_minus_Horizons_00UT_separation_deg":geodesic(obs,15*(h+m/60),de),
                    "Horizons_SOHOCentric_ICRF_RA_Dec": [obs[3],obs[4]],
                    "Horizons_projected_3sigma_arcsec":[float(obs[7]),float(obs[8])],
                    "two_body_L42_vs_N33_00TDB_separation_arcsec":sep})
assert len(results)==4
args.output.write_text(json.dumps({"positions":results,"quicklook_images":len(filenames),"receipts":receipts},
                                  indent=2,allow_nan=False)+"\n")
