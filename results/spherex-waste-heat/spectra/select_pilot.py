"""Independent blind check: pre-selected Gaia pair from the declared blind field."""
import csv
from select_fields import COLUMNS,ROOT,pos,healpix_nest

ids={'4008195525102016000':'pilot','4007948689037185536':'pilot_cmp1'}
rows=[]
for star in csv.DictReader((ROOT.parent/'blind/stars.csv').open()):
    if star['source_id'] in ids:
        star['ref_epoch']='2016.0'
        ra,dec=pos(star)
        rows.append(dict(label=ids[star['source_id']],**{k:star[k] for k in COLUMNS.split(',')},
                         ra_2025_9=ra,dec_2025_9=dec,tile=healpix_nest(5,ra,dec)))
assert len(rows)==len(ids)
with (ROOT/'positions_pilot.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
