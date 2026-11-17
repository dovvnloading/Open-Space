"""Check detector-registered aperture growth and a simultaneous field star at each S110 hit."""
import argparse
import gzip
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from ffi_acquire import propagated, tsv_rows
from ffi_images import event_masks
from latest_ffi import one
from sector110_followup import SOURCE


def image_diagnostics(rows, event):
    t=np.array([r['time'] for r in rows])
    inside, left, right=event_masks(t,event['epoch_btjd'],event['duration_min']/1440.,None)
    outside=left|right
    center=np.median([[r['x'],r['y']] for r in rows if abs(r['time']-event['epoch_btjd'])>.01],axis=0)
    # Stamps follow the target. First register to absolute detector pixels,
    # keeping only pixels available in every event/shoulder frame.
    chosen=[r for i,r in enumerate(rows) if inside[i] or outside[i]]
    origins=np.array([r['stamp_origin'] for r in chosen])
    x0,y0=origins.max(axis=0)
    x1,y1=origins.min(axis=0)+9
    yy,xx=np.mgrid[y0:y1,x0:x1]
    cube=np.array([np.array(r['stamp'])[yy-r['stamp_origin'][1],xx-r['stamp_origin'][0]]-r['sky']
                   for r in chosen])
    intime=np.array([abs(r['time']-event['epoch_btjd'])<=event['duration_min']/2880.+1e-8 for r in chosen])
    shoulders=cube[~intime].mean(axis=0)
    during=cube[intime].mean(axis=0)
    delta=shoulders-during
    dist=np.hypot(xx-center[0],yy-center[1])
    curve=[]
    for radius in (1.5,2.,2.5,3.,3.5):
        ap=dist<=radius
        expected=sum(np.hypot(x-center[0],y-center[1])<=radius for x in range(round(center[0])-4,round(center[0])+5)
                     for y in range(round(center[1])-4,round(center[1])+5))
        if ap.sum()!=expected: continue
        curve.append(dict(radius=radius, shoulder_flux=float(shoulders[ap].sum()),
                          in_flux=float(during[ap].sum()),depth=float(delta[ap].sum()/shoulders[ap].sum())))
    def centroid(image):
        weights=np.where(dist<=3., np.maximum(0.,image), 0.)
        return np.array([(weights*z).sum()/weights.sum() for z in (xx,yy)])
    coords=np.array([[r['x'],r['y']] for r in chosen])
    return dict(inside_frames=[r['filename'] for r,b in zip(chosen,intime) if b],
                flank_frames=[r['filename'] for r,b in zip(chosen,intime) if not b],
                fixed_aperture_center=center.tolist(), centroid_change_pix=(centroid(during)-centroid(shoulders)).tolist(),
                wcs_change_pix=(coords[intime].mean(axis=0)-coords[~intime].mean(axis=0)).tolist(),
                grid_origin_xy=[int(x0),int(y0)], difference_e_s=np.round(delta,3).tolist(), fixed_apertures=curve,
                input_ranges=[{k:r[k] for k in ('filename','start','end','header_sha256','strip_sha256')}
                              for r in chosen])


def main(cache, out):
    series=list(map(json.loads,gzip.open(out/'frames.jsonl.gz','rt')))
    report=json.loads((out/'search.json').read_text())
    neighbours=tsv_rows(cache/f'{SOURCE}-gaia.tsv')
    field=min((r for r in neighbours if r['Source']!=SOURCE),key=lambda r:float(r['Gmag']))
    fc=propagated(field,2026.8)
    fstar=dict(camera=1,ccd=4,ra_obs=fc[0],dec_obs=fc[1])
    fcache=cache/'field_star'
    fcache.mkdir(exist_ok=True)
    def field_frame(name):
        p=fcache/(name+'.json')
        if not p.exists():
            row=one(name,fstar,keep_pixels=True)
            if 'stamp' not in row: raise ValueError(str(row))
            p.write_text(json.dumps(row)+'\n')
        return json.loads(p.read_text())
    results=[]
    for event in report['search']['isolated']:
        selected=[r for r in series if r.get('min_target_fraction',0)>=.8 and abs(r['time']-event['epoch_btjd'])<.03]
        target=[json.loads((cache/'cache'/(r['filename']+'.json')).read_text()) for r in selected]
        with ThreadPoolExecutor(6) as pool:
            comparisons=list(pool.map(field_frame,[r['filename'] for r in selected]))
        results.append(dict(event=event, target=image_diagnostics(target,event),comparison_gaia=field['Source'],
                            comparison_propagated_radec=list(fc),comparison=image_diagnostics(comparisons,event)))
    (out/'image_vet.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cache',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args()
    main(a.cache,a.out)
