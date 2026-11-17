"""List KMT alerts for additional full-curve, baseline-constrained continuous fits."""
import argparse
import csv
import json
from pathlib import Path

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,default=Path.cwd())
    p.add_argument('--threshold',type=float,default=25.);a=p.parse_args()
    root=Path(__file__).resolve().parent
    seen=set()
    for row in csv.DictReader((root/'kmt/short-lead-audit.tsv').open(),delimiter='\t'):
        if row['event'] in seen or row['status']!='ok':continue
        if not row.get('second_site_quality_gain') or float(row['second_site_quality_gain'])<16:continue
        p=a.cache/'staging/bound'/(row['event']+'.json')
        if not p.is_file():continue
        r=json.loads(p.read_text())
        if (row['model']=='ps' and r['ps_te']<.5 and r['ps_gain']>=80) or (
            row['model']=='fs' and (r['fs_rho']>1 or r['fs_te']<.5) and
            r['fs_gain']-max(r['ps_gain'],r['bin_gain'])>a.threshold):
            print(r['event']);seen.add(r['event'])
