"""Fit fixed-sky source positions in one visit, then lock them in the others."""
import argparse,csv
from pathlib import Path
import numpy as np
from wide_survey import visit_inputs,projected,sample,ARCSEC,inverse_sky,sky_xy


def run(files,peaks,out):
    ds=[visit_inputs(p) for p in sorted(files)]
    for d in ds:d['good']=d['good'].astype('f4')
    rows=list(csv.DictReader(open(peaks)))
    results=[]
    for row in rows:
        a=float(row['a_km']);ph=np.deg2rad(float(row['phase_deg']))
        for i,d in enumerate(ds):
            if float(row[f'v{i+1}_score'])<5 or not int(row[f'ok{i+1}']):continue
            o=projected(d,a,np.array([ph]))[0];mid=len(o)//2
            ra,dec=inverse_sky(d['path_ra'][mid],d['path_dec'][mid],o[mid,0]/ARCSEC,o[mid,1]/ARCSEC)
            # The anchor coadd alone determines the refined sidereal position.
            dd=np.arange(-.30,.301,.025);east,north=np.meshgrid(dd,dd)
            ras,decs=inverse_sky(ra,dec,east.ravel()/ARCSEC,north.ravel()/ARCSEC)
            def get(v,r,d):
                off=sky_xy(r[:,None],d[:,None],v['path_ra'],v['path_dec'],1.,(0.,0.));off[:,:,0]*=-1
                return sample(v,off)
            anchor,usable=get(d,ras,decs)
            if not usable.any():continue
            k=np.argmax(np.where(usable,anchor,-np.inf))
            scores=[get(v,ras[k:k+1],decs[k:k+1]) for v in ds]
            results.append([a,row['phase_deg'],str(d['visit']),ras[k],decs[k],
                            *[float(q[0][0]) for q in scores],*[bool(q[1][0]) for q in scores]])
    with open(out,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['a_km','track_phase_deg','fixed_position_fit_visit','ra_icrf_deg','dec_icrf_deg',
                                   'fixed_sky_score_visit1','fixed_sky_score_visit2','fixed_sky_score_visit3','ok1','ok2','ok3'])
        w.writerows(results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('peaks');p.add_argument('out');p.add_argument('maps',nargs='+')
    a=p.parse_args();run(a.maps,a.peaks,a.out)
