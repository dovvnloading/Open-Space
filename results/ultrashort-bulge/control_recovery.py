"""Compare each in-window published identification against the complete common screen."""
import csv
from consolidate import HERE, rows, read_results
from pathlib import Path
import argparse

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,default=Path.cwd());a=p.parse_args()
    catalog={r['known_id']:r for r in rows(HERE/'other/known_roster.csv')}
    stage=a.cache/'staging'
    sets=[('primary grid',read_results(stage/'screen')),('KMT reference-flux grid',read_results(stage/'bound')),
          ('continuous capped free-grid seed',read_results(stage/'crossings'))]
    out=[]
    for alias in rows(HERE/'other/known_matches.csv'):
        if alias['status']!='known by published ID':continue
        n,k=alias['alert_id'],alias['known_id'];pub=catalog[k]
        base=dict(published_id=k,alert_id=n,published_te_d=pub['te_d'],published_theta_uas=pub['theta_e_uas'])
        for method,data in sets:
            if n not in data:continue
            s=data[n];rho=s.get('fs_rho',s.get('rho'))
            out.append(dict(**base,method=method,ps_te_d=s['ps_te'],fs_te_d=s['fs_te'],rho=rho,
                            selected=int(min(s['ps_te'],s['fs_te'])<.5 or rho>1)))
    with (HERE/'tables/known-recovery.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(out[0]),lineterminator='\n');w.writeheader();w.writerows(out)
