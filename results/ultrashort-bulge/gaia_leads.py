"""Audit Gaia time sampling for each catalog alert that meets a numerical grid cut."""
import csv
import json
from pathlib import Path
from consolidate import HERE,rows


def main():
    meta={r['event']:r for r in rows(HERE/'other/gaia.csv')}
    output=[]
    for p in sorted(Path('staging/screen').glob('Gaia*.json')):
        r=json.loads(p.read_text())
        if r['status']!='ok':continue
        if min(r['ps_te'],r['fs_te'])>=.5 and r['fs_rho']<=1:continue
        c=meta[r['event']]
        output.append(dict(event=r['event'],ps_te=r['ps_te'],fs_te=r['fs_te'],fs_rho=r['fs_rho'],
            estimated_gain_fs_minus_ps=round(r['fs_gain']-r['ps_gain'],4),
            n_before_best_single_0p5d=r['n_before_0.5d'],n_after_best_single_0p5d=r['n_after_0.5d'],
            classification=c['classification'],alert_comment=c['comment'],lightcurve_url=c['url']))
    with (HERE/'tables/gaia-flagged.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)


if __name__=='__main__':main()
