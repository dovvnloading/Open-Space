# SWAN observation-time consistency

Official MOST instructions describe matched ephemeris times as the midpoint of an observation: https://irsa.ipac.caltech.edu/applications/MOST/instructions.html .

For all nine retrieved SWAN pilot cutouts, compare `0.5*(mjd_start+mjd_end)` from the frozen MOST manifest against `MJD-AVG` in the FITS IMAGE header. Maximum absolute discrepancy is0.000044seconds, consistent with serialization precision. The table's `start_date` is the exposure start, not the ephemeris epoch; reports use it as an observation identifier only.

This checks archive time consistency, not an independent orbit integration or geocentric-to-spacecraft ephemeris correction. Small December centroid offsets provide a separate descriptive image-domain check, not a substitute for full astrometric modeling.

Reproduce after downloading the pilot with the supplied script:

```python
import json
from examples.spherex_moving_screen import parse_most
rows={r['observationid']+r['productid']:r for r in parse_most('results/spherex-dr1-moving-spectra/most-swan.tbl')}
for x in json.load(open('results/spherex-dr1-moving-spectra/pilot-swan-retry.json')):
    r=rows[x['observationid']+x['detector']]
    print(x['observationid'],(0.5*(float(r['mjd_start'])+float(r['mjd_end']))-x['mjd'])*86400)
```
