"""Match targets to each RACS-low2 detection at its scan epoch and sky-offset controls."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ARCSEC = np.pi / (180 * 3600)
FIELDS = ['sample','target','gaia_dr3_id','dist_pc','astrometry','motion_applied',
          'target_ra_deg','target_dec_deg','radio_id','source_association','source_type',
          'mjd','sbid','separation_arcsec','radius_arcsec','i_peak_mjy_beam',
          'i_rms_mjy_beam','v_peak_mjy_beam','v_rms_mjy_beam','v_over_i_pct',
          'abs_v_over_i_ge20pct','leakage_pct','offset_arcsec','offset_pa_deg']


def frame(ra_deg, dec_deg):
    a,d = np.deg2rad(ra_deg),np.deg2rad(dec_deg)
    central = np.column_stack([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
    east = np.column_stack([-np.sin(a),np.cos(a),np.zeros_like(a)])
    north = np.column_stack([-np.cos(a)*np.sin(d),-np.sin(a)*np.sin(d),np.cos(d)])
    return central,east,north


def write_csv(path, rows):
    with path.open('w',newline='') as stream:
        w=csv.DictWriter(stream, fieldnames=FIELDS, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def match(sample, path, radio):
    s=pd.read_csv(path,dtype={'gaia_dr3_id':str})
    s=s.loc[np.isfinite(s.ra_deg) & np.isfinite(s.dec_deg)].reset_index(drop=True)
    pos,east,north=frame(s.ra_deg,s.dec_deg)
    can_move=s[['ref_epoch','pmra_masyr','pmdec_masyr']].notna().all(axis=1).to_numpy()
    pm=np.where(can_move[:,None],(s.pmra_masyr.to_numpy()[:,None]*east+
                                 s.pmdec_masyr.to_numpy()[:,None]*north)*ARCSEC/1000,0)
    epoch=s.ref_epoch.fillna(2022.4).to_numpy()
    gaia=s.gaia_dr3_id.fillna('').to_numpy()
    # A coarse tree at a central epoch, then exact scan-epoch tangent-plane motion
    # normalised in Cartesian space (mu_alpha* already includes cos(delta)).
    central=pos+(2022.4-epoch)[:,None]*pm
    central/=np.linalg.norm(central,axis=1)[:,None]
    tree=cKDTree(central)
    result=[]; controls=[]
    maxspeed=np.max(np.linalg.norm(pm,axis=1))/ARCSEC
    flags = [(0,0)]+[(r,a) for r in (60,120,300,600) for a in range(0,360,45)]
    for row in radio.itertuples(index=False):
        r,e,n=frame(np.array([row.ra]),np.array([row.dec_corr]))
        yr=2000+(row.scan_start_mjd-51544.5)/365.25
        for dist,angle in flags:
            t=np.deg2rad(angle)
            ref=(np.cos(dist*ARCSEC)*r+np.sin(dist*ARCSEC)*(e*np.sin(t)+n*np.cos(t)))[0]
            window=(row.search_radius+abs(yr-2022.4)*maxspeed+0.05)*ARCSEC
            near=np.array(tree.query_ball_point(ref,window),dtype=int)
            if not len(near): continue
            obs=pos[near]+(yr-epoch[near])[:,None]*pm[near]
            obs/=np.linalg.norm(obs,axis=1)[:,None]
            sep=2*np.arcsin(np.minimum(1,np.linalg.norm(ref-obs,axis=1)/2))/ARCSEC
            for ix in np.flatnonzero(sep <= row.search_radius):
                j=near[ix]; star=s.iloc[j];x,y,z=obs[ix]
                d=[sample,star.hostname,gaia[j],star.dist_pc if pd.notna(star.dist_pc) else '',
                   star.astrometry,bool(can_move[j]),round(np.rad2deg(np.arctan2(y,x))%360,9),
                   round(np.rad2deg(np.arcsin(z)),9),row.id,row.source_association,row.source_type,
                   row.scan_start_mjd,row.sbid,round(sep[ix],4),row.search_radius,row.peak_flux,row.noise,
                   round(row.peak_flux_v,7),row.noise_v,row.peak_vi,abs(row.peak_vi)>=20,row.leakage,dist,angle]
                (controls if dist else result).append(dict(zip(FIELDS,d)))
    return result,controls,dict(sample=sample,rows_with_position=len(s),rows_motion_applied=int(sum(can_move)),
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),target_pairs=len(result),
        distinct_targets=len({x['target'] for x in result}),control_pairs=len(controls),
        expected_chance_pairs=len(controls)/32)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sample_dir',type=Path)
    args=p.parse_args()
    radio=pd.read_csv(HERE/'racs_low2_v.csv')
    results=[];controls=[];summary=[]
    for label in ('hosts','gaia_nearby','cns5','dwarfs'):
        r,c,m=match(label,args.sample_dir/(label+'.csv'),radio)
        results.extend(r);controls.extend(c);summary.append(m)
    write_csv(HERE/'racs_matches.csv',results);write_csv(HERE/'racs_offset_matches.csv',controls)
    (HERE/'racs_control_counts.json').write_text(json.dumps({'radio_records':len(radio),
        'control_positions_per_radio_record':32,'samples':summary},indent=2)+'\n')
