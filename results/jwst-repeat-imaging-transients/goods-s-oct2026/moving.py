import argparse
import urllib.request as u, urllib.parse as p, json
from astropy.coordinates import SkyCoord
import astropy.units as au
from pathlib import Path

parser = argparse.ArgumentParser(description="Two-pass known small-body check with a JWST observer state")
parser.add_argument("--ra", type=float, default=53.214163)
parser.add_argument("--dec", type=float, default=-27.769795)
parser.add_argument("--utc", default="2026-10-10T00:00:00")
parser.add_argument("--output", type=Path, default=Path("known_moving_objects.json"))
args = parser.parse_args()
s = SkyCoord(args.ra, args.dec, unit="deg")
from astropy.time import Time

utc = args.utc
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
    "fov-dec-center": ("M" if s.dec.deg < 0 else "") + abs(s.dec).to_string(
        unit=au.deg, sep="-", precision=1, pad=True
    ),
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
    args.output.open("w"),
    indent=2,
)
