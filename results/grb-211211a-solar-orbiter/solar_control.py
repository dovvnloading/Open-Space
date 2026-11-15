"""Extract the independent GOES-16 X-ray control, verifying the NCEI archive bytes."""
import csv
import hashlib
import sys
from datetime import datetime, timedelta, timezone
from io import BytesIO
from pathlib import Path
from urllib.request import urlopen

import h5py

URL = ('https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/'
       'goes16/l2/data/xrsf-l2-avg1m_science/2021/12/'
       'sci_xrsf-l2-avg1m_g16_d20211211_v2-2-2.nc')
SHA256 = '3a21beea1838dbd74faae36378c7e4f43ba68373c3a631ccb4201b9e30d22191'
START = datetime(2021,12,11,13,0,tzinfo=timezone.utc)
EPOCH = datetime(2000,1,1,12,tzinfo=timezone.utc)


def extract(raw):
    if hashlib.sha256(raw).hexdigest() != SHA256:
        raise ValueError('The GOES-16 input differs from the pinned release')
    with h5py.File(BytesIO(raw)) as data:
        for i, offset in enumerate(data['time'][:]):
            t = EPOCH + timedelta(seconds=float(offset))
            if START <= t < START+timedelta(minutes=26):
                yield [t.isoformat(), *[float(data[key][i]) for key in ('xrsa_flux','xrsb_flux')],
                       *[int(data[key][i]) for key in ('xrsa_flag','xrsb_flag')]]


if __name__ == '__main__':
    if len(sys.argv)>1:
        raw = Path(sys.argv[1]).read_bytes()
    else:
        with urlopen(URL, timeout=60) as source:
            raw = source.read()
    with Path(__file__).with_name('goes16.csv').open('w',newline='') as handle:
        writer = csv.writer(handle,lineterminator='\n')
        writer.writerow(['utc','flux_0p5_4_angstrom_W_m2','flux_1_8_angstrom_W_m2',
                         'short_quality_flag','long_quality_flag'])
        writer.writerows(extract(raw))
