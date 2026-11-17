"""List KMT alerts for additional full-curve, baseline-constrained continuous fits."""
import argparse
import csv
import json
from pathlib import Path


def threshold_crossings(cache):
    """High-signal short free fits crossing half a day after imposing the flux cap."""
    for p in sorted((cache/'staging/bound').glob('*.json')):
        free=cache/'staging/screen'/p.name
        if not free.is_file():continue
        a=json.loads(free.read_text());b=json.loads(p.read_text())
        if a.get('status')!='ok' or b.get('status')!='ok':continue
        if (min(a['ps_te'],a['fs_te'])<.5 and .5<=min(b['ps_te'],b['fs_te'])<=3
                and b['fs_rho']<=1 and a['fs_gain']-max(a['ps_gain'],a['bin_gain'])>40
                and min(a['n_after_0.5d'],a['n_before_0.5d'])>=2):
            yield p.stem


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--threshold',type=float,default=25.)
    p.add_argument('--crossings',action='store_true');a=p.parse_args()
    if a.crossings:
        for event in threshold_crossings(a.cache):print(event)
        raise SystemExit(0)
    root=Path(__file__).resolve().parent
    seen=set()
    for row in csv.DictReader((root/'kmt/short-lead-audit.tsv').open(),delimiter='\t'):
        if row['event'] in seen or row['status']!='ok':continue
        def value(key):return float(row.get(key) or 0)
        independent=value('second_site_quality_gain')>=16
        one_site=(0<value('te_d')<.5 and value('n_before')>=2 and value('n_after')>=2
                  and value('n_gt5')-value('n_gt5_bad_seeing')>=4 and value('quality_gain')>=80
                  and row.get('t0_minus_catalog_d') and abs(value('t0_minus_catalog_d'))<3)
        if not (independent or one_site):continue
        p=a.cache/'staging/bound'/(row['event']+'.json')
        if not p.is_file():continue
        r=json.loads(p.read_text())
        if one_site or (row['model']=='ps' and r['ps_te']<.5 and r['ps_gain']>=80) or (
            row['model']=='fs' and (r['fs_rho']>1 or r['fs_te']<.5) and
            r['fs_gain']-max(r['ps_gain'],r['bin_gain'])>a.threshold):
            print(r['event']);seen.add(r['event'])
