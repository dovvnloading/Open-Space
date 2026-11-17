"""Resolved Gaia neighbors and a local optical source-density control in the same cones."""
import argparse,csv
from pathlib import Path
import numpy as np
from scipy.stats import chi2
def main(out):
    sums=[]
    for host,sid in [('HD53143','5479222240596469632'),('HD38206','2967287596507432704')]:
        rows=list(csv.DictReader((out/(host+'-gaia.csv')).open()));h=next(r for r in rows if r['source_id']==sid)
        r0,d0=float(h['ra']),float(h['dec'])
        # Both 40-arcsec control disks lie fully inside their original 43.2-arcsec queries.
        cq={'HD53143':(104.9963,-61.3344),'HD38206':(85.8404,-18.5576)}[host]
        assert np.hypot((r0-cq[0])*np.cos(np.deg2rad(d0)),d0-cq[1])*3600+40<=43.2
        sep=np.array([3600*np.hypot((float(r['ra'])-r0)*np.cos(np.deg2rad(d0)),float(r['dec'])-d0) for r in rows])
        n=int(((sep>=10)&(sep<40)).sum());area=np.pi*(40**2-10**2);search=np.pi*(5.5**2-2.5**2)
        sums.append(dict(host=host,Catalogue='Gaia DR3 all G',control_inner_arcsec=10,control_outer_arcsec=40,count=n,control_area_sq_arcsec=area,
                         search_annulus_expected=n*search/area,poisson_lower95=.5*chi2.ppf(.025,2*n)*search/area if n else 0,
                         poisson_upper95=.5*chi2.ppf(.975,2*n+2)*search/area,observed_search_annulus=int(((sep>=2.5)&(sep<5.5)).sum())))
    with (out/'gaia_control_counts.csv').open('w') as f:w=csv.DictWriter(f,list(sums[0]));w.writeheader();w.writerows(sums)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);main(p.parse_args().out)
