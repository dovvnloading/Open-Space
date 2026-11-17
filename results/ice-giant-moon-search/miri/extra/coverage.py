"""Map propagated circular-orbit samples into all 126 independent native MIRI footprints."""
import argparse,csv
from pathlib import Path
import numpy as np
from astropy.time import Time
from astropy.wcs import WCS
from fetch import ROOT,products
from geometry import PLANET,AU,C,ARCSEC,basis,frequency,tangent_basis,direction
from miri_audit import headers,write_csv

OUT=Path(__file__).resolve().parent

def get_states():
    rows=list(csv.DictReader((OUT/'positions.csv').open()))
    lookup={}
    for r in rows: lookup.setdefault(r['file'],{})[int(r['body'])]=np.array([float(r[x]) for x in ['ra_deg','dec_deg','diameter_arcsec','r_au','delta_au']])
    return lookup


def emission(f,state):
    t=Time(float(f['exp_mid_mjd_utc']),format='mjd',scale='utc').tdb.jd
    return t-state[4]*AU/C/86400


def frame(f,product,cache,states,receipts):
    state=states[PLANET[f['planet']]]
    _,ph,h,_=headers(product,cache,receipts);w=WCS(h,naxis=2);nx=h['NAXIS1'];ny=h['NAXIS2']
    los,east,north=tangent_basis(*state[:2]);u,v,_=basis(f['planet'],emission(f,state))
    # Broad outer box only for selecting WCS inverse points. Final test uses SIP pixels.
    edges=np.linspace(0,1,65);x=np.r_[edges*(nx-1),edges*0+nx-1,(1-edges)*(nx-1),edges*0];y=np.r_[edges*0,edges*(ny-1),edges*0+ny-1,(1-edges)*(ny-1)]
    ra,dec=w.all_pix2world(x,y,0)
    dra=(ra-state[0]+180)%360-180
    return dict(f=f,w=w,nx=nx,ny=ny,state=state,states=states,pixar_sr=h['PIXAR_SR'],t=emission(f,state),u=u,v=v,
                dist=state[4]*AU,los=los,skybox=(dra.min()-.0001,dra.max()+.0001,dec.min()-.0001,dec.max()+.0001))


def pixels(g,radii,phase,origin):
    theta=phase+frequency(g['f']['planet'],radii)*(g['t']-origin)*86400
    xyz=np.cos(theta)[:,None]*g['u']+np.sin(theta)[:,None]*g['v']
    line=g['dist']*g['los']+radii[:,None]*xyz
    ra=np.rad2deg(np.arctan2(line[:,1],line[:,0]))%360
    dec=np.rad2deg(np.arctan2(line[:,2],np.hypot(line[:,0],line[:,1])))
    dra=(ra-g['state'][0]+180)%360-180;l,r,b,t=g['skybox']
    keep=(dra>=l)&(dra<=r)&(dec>=b)&(dec<=t);indices=np.flatnonzero(keep)
    x,y=g['w'].all_world2pix(ra[keep],dec[keep],0,quiet=True) if len(indices) else (np.array([]),np.array([]))
    eligible=np.isfinite(x+y)&(x>=0)&(x<=g['nx']-1)&(y>=0)&(y<=g['ny']-1)
    return indices[eligible], x[eligible], y[eligible]


def run(cache):
    products_ = products(); states=get_states();bygroup={};records=[];r_rows=[];visits=[]
    obs={r['obsid']:r for r in csv.DictReader((ROOT/'manifests/observations.csv').open())}
    for planet,outer in [('Uranus',650000),('Neptune',300000)]:
        selected=[(f,p) for f,p in products_ if f['planet']==planet];frames=[]
        for f,p in selected:
            g=frame(f,p,cache,states[f['file']],{});g['group']=[obs[k]['obs_id'] for k in p['parent_obsids'].split('|') if obs[k]['instrument_name']=='MIRI/IMAGE'][0];frames.append(g)
        origin=min(g['t'] for g in frames);rs=np.arange(50000,outer+1,5000);phase=np.deg2rad(np.arange(0,360,1.))
        rr=np.repeat(rs,len(phase));ph=np.tile(phase,len(rs));checks={}
        for g in frames:
            idx,_,_=pixels(g,rr,ph,origin);vis=g['f']['visit'];k=(vis,g['group'])
            if k not in checks:checks[k]=np.zeros(len(rr),int)
            checks[k][idx]+=1
            records.append(dict(file=g['f']['file'],group=g['group'],planet=planet,visit=vis,filter=g['f']['filter'],planet_on_chip=g['f']['planet_on_chip'],
                                planet_to_chip_arcsec=g['f']['planet_to_chip_arcsec'],sampled_cells=len(rr),geometric_cells=len(idx),
                                minimum_sampled_km=np.min(rr[idx]) if len(idx) else '',maximum_sampled_km=np.max(rr[idx]) if len(idx) else ''))
        for (v,group),counts in checks.items():
            use=[g for g in frames if g['group']==group]
            for j,r in enumerate(rs):
                c=counts.reshape(-1,len(phase))[j]
                r_rows.append(dict(planet=planet,group=group,visit=v,a_km=r,n_phase=len(c),at_least_one=int((c>0).sum()),at_least_two=int((c>=2).sum())))
            ec=counts.reshape(-1,len(phase)).any(axis=1)
            bygroup[group]=dict(group=group,planet=planet,visit=v,filter=use[0]['f']['filter'],exposures=len(use),planet_centers_on_chip=sum(int(g['f']['planet_on_chip']) for g in use),
                                smallest_sampled_radius_km=rs[ec].min() if ec.any() else '',largest_sampled_radius_km=rs[ec].max() if ec.any() else '',
                                sampled_orbit_cells=len(rr),at_least_two_exposures=int((counts>=2).sum()),epoch_jd_tdb=origin)
        for i,(key,ct) in enumerate(checks.items()):
            for (other,oc) in list(checks.items())[i+1:]:
                if key[0]==other[0]:continue
                valid=(ct>=2)&(oc>=2);visits.append(dict(planet=planet,group_1=key[1],group_2=other[1],shared_orbit_cells=int(valid.sum()),
                    min_radius_km=rr[valid].min() if valid.any() else '',max_radius_km=rr[valid].max() if valid.any() else '',epoch_jd_tdb=origin))
    for name,rows in [('exposure_coverage',records),('group_coverage',list(bygroup.values())),('radius_coverage',r_rows),('visit_overlaps',visits)]:
        write_csv(OUT/(name+'.csv'),list(rows[0]),rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('header_cache',type=Path);run(p.parse_args().header_cache)
