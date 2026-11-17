"""Finite host-distance colour test and empirical foreground/background star counts."""
import csv,hashlib,io,json,sys,tarfile
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
def main(path):
    receipt=next(r for r in json.loads((HERE.parent/'miri/calibration_inputs.json').read_text()) if r['filename']=='bobcat.tar.gz')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=receipt['sha256']:raise ValueError(path)
    distance=1000/175.3131
    # Model log10(mJy at 10pc). Do not fit a free flux scale at known distance.
    pred=[]
    with tarfile.open(path) as t:
        for z in ('-0.5','+0.0','+0.5'):
            data=t.extractfile('photometry_tables/flux_table_JWST'+z).read().decode().splitlines()
            terms=next(s.split()[8:] for s in data if s.lstrip().startswith('Teff'))
            for line in data:
                if '*' in line:continue
                cols=line.split()
                try: temp=float(cols[0]);g=float(cols[1])
                except (ValueError,IndexError):continue
                f10,f15=[1e3*10**float(cols[6+terms.index(b)])*(10/distance)**2 for b in ('F1065C','F1550C')]
                if abs(f10-89.85)<=3*8.55:
                    pred.append(dict(metallicity=z,teff=temp,logg=g,F1065C_uJy=f10,F1550C_uJy=f15))
    with (HERE/'field_model_grid.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(pred[0]));w.writeheader();w.writerows(pred)
    with (HERE.parent/'stellar/hip25878/gaia.csv').open() as f:g=list(csv.DictReader(f))
    n=sum(1 for r in g if r['source_id']!='3209938366665770752' and float(r['phot_g_mean_mag'])<=17.)
    result=dict(host_distance_pc=distance,f1065c_uJy=89.85,ring_noise_uJy=8.55,
                f1550c_forced_uJy=30.53,f1550c_ring_noise_uJy=15.13,
                qualifying_models=len(pred),f1550c_model_range_uJy=[min(v['F1550C_uJy'] for v in pred),max(v['F1550C_uJy'] for v in pred)],
                gaia_cone_radius_arcsec=.02222*3600,G_limit=17.,gaia_count_excluding_host=n,
                gaia_density_per_arcmin2=n/(np.pi*(.02222*60)**2))
    (HERE/'field_color.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main(Path(sys.argv[1]))
