"""Summarize observed sampling in the 27 KMT fields, using the parent cadence pipeline."""
import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
import numpy as np
from freeze import HERE


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cadence',type=Path,default=HERE.parent/'tables/cadence.csv.gz')
    a=p.parse_args()
    opener=gzip.open if a.cadence.suffix=='.gz' else open
    with opener(a.cadence,'rt') as h:
        coverage={r['event']:r for r in csv.DictReader(h) if not r['error'] and r['event'].startswith('KMT-')}
    table=defaultdict(list)
    for year in (2024,2025):
        with (HERE/f'catalog-{year}.tsv').open() as h:
            for r in csv.DictReader(h,delimiter='\t'):
                field=re.search(r'\d{2}',r['Field']).group()
                table[(year,field)].append(coverage.get(r['Event']))
    rows=[]
    for (year,field),events in sorted(table.items()):
        data=[r for r in events if r]
        rec=dict(year=year,field=field,listed=len(events),with_measured_sampling=len(data))
        if data:
            rec['median_exposures_per_event']=round(float(np.median([float(r['distinct_epochs_in_season']) for r in data])),1)
            for s in ['0.03','0.1','0.25','0.5','1.5']:
                fractions=[float(r['days_'+s])/float(r['season_d']) for r in data if float(r['season_d'])>0]
                if fractions:
                    rec[f'median_fraction_plusminus_{s}d']=round(float(np.median(fractions)),5)
        rows.append(rec)
    with (HERE/'field-coverage.tsv').open('w') as h:
        w=csv.DictWriter(h,list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)


if __name__=='__main__':main()
