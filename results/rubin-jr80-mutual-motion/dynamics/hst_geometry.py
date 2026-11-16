"""Fetch the HST-epoch sightline and range of (556068) from JPL Horizons."""

from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
params = dict(
    format="text",
    COMMAND="'556068;'",
    CENTER="'@-48'",
    EPHEM_TYPE="'OBSERVER'",
    TLIST="'2461205.091'",
    QUANTITIES="'1,20'",
    ANG_FORMAT="'DEG'",
    CSV_FORMAT="'YES'",
    EXTRA_PREC="'YES'",
    OBJ_DATA="'NO'",
)
url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
if __name__ == "__main__":
    with urlopen(url, timeout=60) as response:
        body = response.read()
    if b"$$SOE" not in body or b"$$EOE" not in body:
        raise RuntimeError(body.decode()[:2000])
    (ROOT / "horizons-hst.txt").write_bytes(body)
