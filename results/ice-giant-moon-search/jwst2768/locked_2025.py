"""Predict Jan 2024 with each pre-existing 2025 circular alias, without refitting it."""
import argparse,csv,json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import generate
from search import Cuts,load_npz,prepare,sample_paths,score,validate_cube


def run(paths,aliases,out):
    out=Path(out);prior=list(csv.DictReader((Path(aliases)/'circular_aliases.csv').open()))
    meta=json.loads((Path(aliases)/'summary.json').read_text());epoch=meta['epoch2']
    par=np.zeros((len(prior),6));par[:,0]=[float(q['a_km']) for q in prior];par[:,-1]=[float(q['phase_late_deg']) for q in prior]
    recs=[dict(cycle=int(q['cycle']),a_km=q['a_km'],phase_deg=q['phase_late_deg'],
                training_feb_rms_px=q['feb_rms'],training_oct_rms_px=q['oct_rms']) for q in prior]
    residuals=[]
    with TemporaryDirectory() as tmp:
        elements=Path(tmp)/'elements.csv';np.savetxt(elements,par,delimiter=',',comments='',header='a_km,e,i_deg,node_deg,peri_deg,mean_deg')
        for fname in paths:
            c=load_npz(fname);band=str(c['band']).split('/')[1];tr=Path(tmp)/'tracks.npz'
            generate(fname,tr,elements=elements,model='circular',epoch=epoch);t=load_npz(tr)
            cuts=Cuts(psf_support=.8);maps=prepare(c,validate_cube(c),cuts)
            n,w=sample_paths(maps,t['x_path'],t['y_path'],cuts);fit=score(n,w,c['visit'],cuts)
            ast=list(csv.DictReader((out/f'{band}_native.csv').open()))
            ix=[list(c['exposure_id']).index(a['frame'].rsplit('_',2)[0]) for a in ast]
            xs=np.array([float(a['sky_x']) for a in ast]);ys=np.array([float(a['sky_y']) for a in ast])
            sep=np.hypot(xs-t['x'][:,ix],ys-t['y'][:,ix]);residuals.append(sep)
            for i,r in enumerate(recs):
                r[band+'_native_rms_px']=float(np.mean(sep[i]**2)**.5)
                r[band+'_native_median_px']=float(np.median(sep[i]))
                r[band+'_pooled_score']=float(fit['snr'][i]);r[band+'_rank']=float(fit['rank'][i])
                for j,v in enumerate(np.unique(c['visit'])):
                    r[band+'_'+v+'_score']=float(fit['visit_snr'][i,j]);r[band+'_'+v+'_n']=int(fit['n'][i,j])
    rms=np.sqrt(np.mean(np.concatenate(residuals,axis=1)**2,axis=1))
    for i,r in enumerate(recs):r['all_2024_native_rms_px']=float(rms[i])
    with (out/'locked_2025_aliases.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=recs[0]);w.writeheader();w.writerows(recs)
    print('2025 top',recs[0]['cycle'],rms[0],'closest to 2024',recs[rms.argmin()]['cycle'],rms.min())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('aliases');p.add_argument('output');p.add_argument('cubes',nargs='+')
    a=p.parse_args();run(a.cubes,a.aliases,a.output)
