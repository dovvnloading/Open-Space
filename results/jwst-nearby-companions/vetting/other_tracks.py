"""Test Gaia and dated Pan-STARRS positions against the MIRI detections."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from archival_images import ensure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FIELDS={
    'hip':dict(host='3209938366665770752',candidate='3209938199163908352',source=(82.87071109,-3.69122188),
               epoch='2025-09-21T22:45:42',relative=(4.019,3.661),ps1='103570828707811060'),
    'g29':dict(host='2660358032257156736',candidate=None,source=(352.193797,5.246759),
               epoch='2025-11-13T19:06:28',relative=(-6.466,1.015),ps1='114293521937476539'),
}

def parpos(s,t):
    ra,dec=float(s['ra']),float(s['dec'])
    a,d=np.deg2rad([ra,dec]);earth=get_body_barycentric('earth',t).xyz.to_value('au')
    tangent=np.array([[-np.sin(a),np.cos(a),0],[-np.cos(a)*np.sin(d),-np.sin(a)*np.sin(d),np.cos(d)]])
    pm=np.array([float(s['pmra']),float(s['pmdec'])]);shift=(t.jyear-2016)*pm-float(s['parallax'])*(tangent@earth)
    return np.array([ra+shift[0]/(3.6e6*np.cos(d)),dec+shift[1]/3.6e6])

def offset(q,p):
    return (np.array(q)-np.array(p))*np.array([np.cos(np.deg2rad(p[1])),1.])*3600

def main():
    par=argparse.ArgumentParser(description=__doc__);par.add_argument('cache',type=Path);arg=par.parse_args()
    for row in json.loads((HERE/'other_archive_inputs.json').read_text()):
        ensure(row,arg.cache)
    results={};epochs=[]
    for tag,setup in FIELDS.items():
        base=arg.cache/(tag+'-archive')
        with (base/'gaia.csv').open() as f:stars={s['source_id']:s for s in csv.DictReader(f)}
        host=stars[setup['host']];date=Time(setup['epoch'])
        hp=parpos(host,date);source=setup['source'];rel=np.array(setup['relative'])
        with (base/'ps1.csv').open() as f:ps=next(r for r in csv.DictReader(f) if r['objID']==setup['ps1'])
        t=fits.getdata(base/'tractor.fits',1);sep=np.hypot((t['ra']-source[0])*np.cos(np.deg2rad(source[1])),t['dec']-source[1])*3600
        i=int(np.argmin(sep))
        result=dict(sky_deg=list(source),relative_arcsec=rel.tolist(),ps1_object=ps['objID'],
            ps1_id_to_miri_arcsec=float(np.linalg.norm(offset([float(ps['raMean']),float(ps['decMean'])],source))),
            legacy_dr10_id=f"{t['brickname'][i]}-{t['objid'][i]}",legacy_to_miri_arcsec=float(sep[i]),
            legacy_morphology=t['type'][i],gaia_host=host)
        if setup['candidate']:
            c=stars[setup['candidate']];gp=parpos(c,date);pred=offset(gp,hp)
            result.update(gaia_field_star=c,gaia_relative_pred_arcsec=pred.tolist(),
                          gaia_relative_resid_arcsec=float(np.linalg.norm(pred-rel)),
                          gaia_absolute_resid_arcsec=float(np.linalg.norm(offset(gp,source))))
        with (base/'ps_epochs.csv').open() as f:dets=list(csv.DictReader(f))
        for dr in sorted(dets,key=lambda r:float(r['obsTime'])):
            dt=Time(float(dr['obsTime']),format='mjd');shift=offset(parpos(host,dt),hp)
            observed=offset([float(dr['ra']),float(dr['dec'])],source)
            epochs.append(dict(field=tag,objID=dr['objID'],detectID=dr['detectID'],date=dt.iso[:10],mjd=dt.mjd,
                ra=dr['ra'],dec=dr['dec'],psfQf=dr['psfQf'],fixed_resid_arcsec=np.linalg.norm(observed),
                host_motion_resid_arcsec=np.linalg.norm(observed-shift),
                host_prediction_shift_arcsec=np.linalg.norm(shift)))
        results[tag]=result
    (HERE/'field_identifications.json').write_text(json.dumps(results,indent=2)+'\n')
    with (HERE/'other_ps1_epochs.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(epochs[0]));w.writeheader();w.writerows(epochs)
if __name__=='__main__':main()
