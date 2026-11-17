"""Sensitivity to Gaia's GSP-Phot distance posterior for the two available candidates."""
from pathlib import Path
import csv
import sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent.parent))
from screen_desi import cartesian
def rows(path):
    with path.open() as f: return list(csv.DictReader(f))
gaia={r['source_id']:r for r in rows(ROOT.parent.parent/'gaia_queue.csv')+rows(ROOT/'gaia.csv')}
gsp=rows(ROOT/'gspphot.csv')
rv={'1717410430709663616':501.7387097558274,'2752545725694798336':-95.19583506584007}
with (ROOT/'gspphot_speeds.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['source_id','gsp_phot_pc','rv_kms','nominal_vgal_kms'])
    for r in gsp:
        sid=r['source_id']
        if sid not in rv: continue
        a=gaia[sid]; dist=float(r['distance_gspphot'])
        args=[a['ra'],a['dec'],dist/1000,a['pmra'],a['pmdec'],rv[sid]]
        _,v=cartesian(*[np.array([float(x)]) for x in args])
        w.writerow([sid,dist,rv[sid],np.linalg.norm(v)])
