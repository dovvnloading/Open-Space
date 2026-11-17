"""Extract the two DESI distance VACs and published quality flags for four inspected spectra."""
from pathlib import Path
import csv,sys
from astropy.io import fits
import numpy as np
ids=[1717410430709663616,4444095982221534720,2752545725694798336,6340918677332110720]
root=Path(__file__).resolve().parent
with fits.open(sys.argv[1],memmap=True) as hf:
    data=hf[1].data; keep=np.isin(data['SOURCE_ID'],ids)
    with (root/'specdis_details.csv').open('w') as f:
        wr=csv.writer(f);wr.writerow(data.names)
        for row in data[keep]:wr.writerow(list(row))
with fits.open(sys.argv[2],memmap=True) as hf:
    data=hf[1].data;keep=np.isin(data['SOURCE_ID'],ids)
    with (root/'spdist_details.csv').open('w') as f:
        wr=csv.writer(f);wr.writerow(data.names+['dist16_kpc','dist50_kpc','dist84_kpc'])
        for row in data[keep]:
            d=[10**((row['G0']-row[k]-10)/5) for k in ('MG_84','MG_50','MG_16')]
            wr.writerow(list(row)+d)
            print(str(row['SOURCE_ID']),*[round(v,3) for v in d],row['FLAG_GOOD'],row['QFLAG_INPUT'])
