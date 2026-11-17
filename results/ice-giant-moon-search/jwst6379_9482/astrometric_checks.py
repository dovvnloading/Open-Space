"""Test a direct SPK calculation, ephemeris longitude offsets and orbit-fit sensitivity."""
import argparse,csv,json
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import least_squares
from astropy.time import Time
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base));sys.path.insert(0,str(base/'jwst9482'))
from link_visits import geometry,search
from displaced import model
from moon_spk import offsets


def visit(cube_path,table,spk):
    c=search.load_npz(cube_path)
    r=list(csv.DictReader(table.open()))
    t=np.array([float(v['mjd_utc']) for v in r])
    xy=np.array([[float(v['x_sky']),float(v['y_sky'])] for v in r])
    ids=np.array([v['frame'] for v in r]);m=np.array([np.flatnonzero(c['exposure_id']==v).item() for v in ids])
    k=np.flatnonzero(c['known_id']==725).item()
    def interp(field):
        return np.array([np.interp(t[i],c['path_mjd'][j],field[j]) for i,j in enumerate(m)])
    ra,dec,delta=[interp(c['path_'+v]) for v in ['ra','dec','delta']]
    known=np.column_stack([interp(c['known_path'][k,:,:,a]) for a in range(2)])
    p=np.column_stack([ra,dec,np.zeros(len(t)),np.zeros(len(t)),delta])
    xyz=offsets(spk,Time(t,format='mjd',scale='utc').tdb.jd,p,target=725)
    sr=np.rad2deg(np.arctan2(xyz[:,1],xyz[:,0]))%360
    sd=np.rad2deg(np.arcsin(xyz[:,2]/np.linalg.norm(xyz,axis=1)))
    spk_xy=geometry.sky_xy(sr,sd,ra,dec,float(c['scale_arcsec']),c['center_xy'])
    ephemcheck=np.linalg.norm(known-spk_xy,axis=1)
    los,e,n=geometry.tangent_basis(ra,dec);u,v,h=geometry.basis('Uranus',2460955.765061979)
    def lift(coords):
        ray=los-(coords[:,0]-c['center_xy'][0])[:,None]*e*c['scale_arcsec']/geometry.ARCSEC+(coords[:,1]-c['center_xy'][1])[:,None]*n*c['scale_arcsec']/geometry.ARCSEC
        arr=delta[:,None]*geometry.AU*((los@h/(ray@h))[:,None]*ray-los)
        return np.stack([arr@u,arr@v],axis=1)
    s=lift(xy);ref=lift(known)
    angles=np.angle((s[:,0]+1j*s[:,1])/(ref[:,0]+1j*ref[:,1]),deg=True)
    rad=np.linalg.norm(s,axis=1)-np.linalg.norm(ref,axis=1)
    dt=(Time(t,format='mjd',scale='utc').tdb.jd-delta*geometry.AU/geometry.C/86400-2460955.765061979)*86400
    def pred(par):return model(par[0],par[1],'Uranus',dt,c,(u,v,h),los,e,n,delta)
    stats=dict(n=len(xy),SPK=search.file_hash(spk),spk_vs_Horizons_median_px=float(np.median(ephemcheck)),
               spk_vs_Horizons_max_px=float(ephemcheck.max()),phase_diff_deg=float(np.median(angles)),
               phase_sigma_deg=float(np.std(angles)),radius_difference_median_km=float(np.median(rad)),
               radius_diff_mad_km=float(np.median(abs(rad-np.median(rad)))),
               nominal_center_offset_median_px=float(np.median(np.linalg.norm(xy-known,axis=1))))
    return xy,ids,pred,stats


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['first_cube','late_cube','first_ast','late_ast','first_spk','late_spk','output']:p.add_argument(n,type=Path)
    a=p.parse_args()
    visits=[visit(a.first_cube,a.first_ast,a.first_spk),visit(a.late_cube,a.late_ast,a.late_spk)]
    groups=np.concatenate([v[1] for v in visits]);allxy=np.concatenate([v[0] for v in visits])
    def residual(par):return np.concatenate([v[2](par) for v in visits])-allxy
    x0=[76418.12,np.deg2rad(254.7467)]
    def fit(mask):
        return least_squares(lambda q:residual(q)[mask].ravel(),x0,x_scale=[.2,.005],loss='soft_l1',f_scale=.6).x
    full=fit(np.ones(len(groups),bool));jk=np.stack([fit(groups!=g) for g in np.unique(groups)])
    sig=np.sqrt((len(jk)-1)/len(jk)*np.sum((jk-jk.mean(0))**2,axis=0))
    report=dict(first=visits[0][3],late=visits[1][3],a_km=full[0],phase_late_deg=np.rad2deg(full[1]),
                n_exposure_blocks=len(jk),conditional_jackknife_sigma_a_km=sig[0],conditional_jackknife_sigma_phase_deg=np.rad2deg(sig[1]),
                leave_one_block_radius_range_km=[jk[:,0].min(),jk[:,0].max()],
                model='The 387-cycle circular equatorial model fixes pole, J2/J4, centering and cycle count. Other cycle solutions are tabulated separately.')
    a.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
