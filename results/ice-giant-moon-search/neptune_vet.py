"""Compare first-pass Neptune tracks to fixed positions in the Proteus and sidereal frames."""
import argparse
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from geometry import ARCSEC, generate, sky_xy
from neptune_scan import photometry
from prepare import inverse_sky
from search import Cuts, load_npz, prepare, sample_paths, validate_cube


def run(paths,survey,out):
    survey=Path(survey);r=list(csv.DictReader((survey/'passing.csv').open()))
    if not r:
        return
    epoch=json.loads((survey/'parameters.json').read_text())['epoch_jd_tdb']
    par=np.zeros((len(r),6));par[:,0]=[float(v['a_km']) for v in r];par[:,-1]=[float(v['phase_deg']) for v in r]
    target_rows=[dict(a_km=v['a_km'],phase_deg=v['phase_deg']) for v in r]; cut=Cuts(psf_support=.8)
    with TemporaryDirectory() as tmp:
        el=Path(tmp)/'elements.csv';np.savetxt(el,par,delimiter=',',header='a_km,e,i_deg,node_deg,peri_deg,mean_deg',comments='')
        datasets=[]
        for path in paths:
            c=load_npz(path);fn=Path(tmp)/'bank.npz';generate(path,fn,elements=el,model='circular',epoch=epoch)
            t=load_npz(fn);p=c['known_path'][list(c['known_id']).index(808)]
            datasets.append((c,t,p))
        # Single constant offset and single fixed ICRF position per hypothesis, shared by both filters.
        dx=np.median(np.concatenate([t['x_path']-p[None,:,:,0] for _,t,p in datasets],axis=1),axis=(1,2))
        dy=np.median(np.concatenate([t['y_path']-p[None,:,:,1] for _,t,p in datasets],axis=1),axis=(1,2))
        stars=[]
        for c,t,p in datasets:
            e=-(t['x_path']-c['center_xy'][0])*float(c['scale_arcsec'])/ARCSEC
            nu=(t['y_path']-c['center_xy'][1])*float(c['scale_arcsec'])/ARCSEC
            # Relative gnomonic coordinates per hypothesis / frame / time.
            ras=[];decs=[]
            for q in range(len(e)):
                ra,dec=inverse_sky(c['path_ra'],c['path_dec'],e[q],nu[q]);ras.append(ra);decs.append(dec)
            stars.append((np.array(ras),np.array(decs)))
        sr,sd=(np.median(np.concatenate([rr[k] for rr in stars],axis=1),axis=(1,2)) for k in [0,1])
        for c,t,p in datasets:
            band=str(c['band']).split('/')[1]
            maps=prepare(c,validate_cube(c,min_visits=1),cut)
            dist=np.hypot(t['x_path']-p[None,:,:,0],t['y_path']-p[None,:,:,1])
            dxs=p[None,:,:,0]+dx[:,None,None];dys=p[None,:,:,1]+dy[:,None,None]
            n,w=sample_paths(maps,dxs,dys,cut);count,s,_,_=photometry(n,w)
            sid=sky_xy(sr[:,None,None],sd[:,None,None],c['path_ra'],c['path_dec'],float(c['scale_arcsec']),c['center_xy'])
            n2,w2=sample_paths(maps,sid[:,:,:,0],sid[:,:,:,1],cut);cnt2,ss,_,_=photometry(n2,w2)
            for i,v in enumerate(target_rows):
                for name,val in [('min_separation_px',dist[i].min()),('median_separation_px',np.median(dist[i])),
                       ('max_separation_px',dist[i].max()),('fixed_proteus_offset_score',s[i]),('fixed_proteus_n',count[i]),
                       ('fixed_sidereal_score',ss[i]),('fixed_sidereal_n',cnt2[i])]:
                    v[band+'_'+name]=val
    with open(out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=target_rows[0]);w.writeheader();w.writerows(target_rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['f140','f210','screen','output']:
        p.add_argument(arg)
    a=p.parse_args();run([a.f140,a.f210],a.screen,a.output)
