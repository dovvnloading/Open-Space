import urllib.request as u, urllib.parse as p, json, re
from astropy.coordinates import SkyCoord
import astropy.units as au
from pathlib import Path

s = SkyCoord(53.214163, -27.769795, unit="deg")
jd = "2461323.5"  # overwritten by astropy Time in UTC/TD below
from astropy.time import Time

utc = "2026-10-10T00:00:00"
jd = str(Time(utc, scale="utc").tdb.jd)
params = {
    "format": "json",
    "COMMAND": "'-170'",
    "CENTER": "'500@399'",
    "EPHEM_TYPE": "'VECTORS'",
    "TLIST": f"'{jd}'",
    "OUT_UNITS": "'KM-S'",
    "REF_PLANE": "'FRAME'",
    "REF_SYSTEM": "'ICRF'",
    "VEC_TABLE": "'2'",
    "CSV_FORMAT": "'YES'",
    "OBJ_DATA": "'NO'",
}
hor = "https://ssd.jpl.nasa.gov/api/horizons.api?" + p.urlencode(params)
r = json.load(u.urlopen(hor, timeout=70))
txt = r.get("result", "")
vec = txt.split("$$SOE")[1].split("$$EOE")[0].strip().split(",")[2:8]
vec = ",".join(t.strip() for t in vec)
q = {
    "xobs": vec,
    "obs-time": utc,
    "fov-ra-center": s.ra.to_string(unit=au.hour, sep="-", precision=2, pad=True),
    "fov-dec-center": s.dec.to_string(
        unit=au.deg, sep="-", precision=1, pad=True
    ).replace("-", "M", 1),
    "fov-ra-hwidth": ".0166667",
    "fov-dec-hwidth": ".0166667",
    "mag-required": "false",
    "two-pass": "true",
    "suppress-first-pass": "true",
    "req-elem": "false",
}
url = "https://ssd-api.jpl.nasa.gov/sb_ident.api?" + p.urlencode(q)
print(url, flush=True)
res = json.load(u.urlopen(url, timeout=240))
print(
    {k: v for k, v in res.items() if k not in ["fields_second", "summary"]}, flush=True
)
json.dump(
    {"horizons_url": hor, "query": url, "horizons": r, "identification": res},
    open("known_moving_objects.json", "w"),
    indent=2,
)
