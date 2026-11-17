"""Compare MER major axes to bright Gaia DR3 stars and 90-degree rotated axes."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent

def unit(ra,dec):
    a,d=np.radians(ra),np.radians(dec)
    return np.c_[np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)]

def transform(stars,coords,ra,dec,theta):
    near=stars.query_ball_point(unit([ra],[dec])[0],2*np.sin(np.radians(120/3600)/2))
    if not near:return {'nearby_stars':0,'aligned':False,'rotated':False,'coincident':False}
    s=coords[near]
    # Local gnomonic offsets are sufficient on 120 arcsec fields.
    de=np.radians(s[:,1]);d=np.radians(dec);a=np.radians(s[:,0]-ra)
    den=np.sin(d)*np.sin(de)+np.cos(d)*np.cos(de)*np.cos(a)
    e=np.rad2deg(np.cos(de)*np.sin(a)/den)*3600
    n=np.rad2deg((np.cos(d)*np.sin(de)-np.sin(d)*np.cos(de)*np.cos(a))/den)*3600
    ang=np.radians(theta)
    # Q1 north-up, east-left mosaic grid. Account for tangent-point rotation separately below.
    ax=-e*np.cos(ang)+n*np.sin(ang);perp=e*np.sin(ang)+n*np.cos(ang)
    distance=np.hypot(e,n);off=distance>=3
    ok=off&(np.abs(perp)<=1.0);rot=off&(np.abs(ax)<=1.0)
    ans=dict(nearby_stars=len(near),coincident=bool((distance<=1).any()),aligned=bool(ok.any()),rotated=bool(rot.any()))
    for selection,label,cross in [(ok,'axis',perp),(rot,'rotated_axis',ax)]:
        if selection.any():
            i=np.where(selection)[0][np.argmin(abs(cross[selection]))]
            ans.update({label+'_gaia_row':int(near[i]),label+'_star_g':float(s[i,2]),label+'_sep_arcsec':float(distance[i]),
                        label+'_cross_arcsec':float(abs(cross[i]))})
    return ans

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stars',type=Path,required=True);p.add_argument('--selected',type=Path,default=HERE/'hats_selected.csv')
    a=p.parse_args();rows=list(csv.DictReader(a.selected.open()))
    tiles={r['tile']:r for r in csv.DictReader((HERE.parent/'euclid/q1_tiles.csv').open())}
    out=[]
    for f in sorted({x['field'] for x in tiles.values()}):
        with gzip.open(a.stars/('gaia_'+f+'.csv.gz'),'rt') as stream:gaia=list(csv.DictReader(stream))
        coords=np.array([[float(s[k]) for k in ('ra','dec','phot_g_mean_mag')] for s in gaia]);stars=cKDTree(unit(coords[:,0],coords[:,1]))
        for r in rows:
            tile=tiles[str(r['tileid'])]
            if tile['field']!=f:continue
            ra,dec,theta=map(float,[r['ra'],r['dec'],r['mer_position_angle']])
            # At RA east of the reference meridian, projected grid north tilts east.
            # Evaluate the projected position angle from the exact inverse TAN Jacobian.
            from math import sin,cos,radians,degrees,atan2
            d0=radians(float(tile['dec_deg']));d=radians(dec);da=radians(ra-float(tile['ra_deg']))
            # tan-plane x=east, y=north; source axes in east-left FITS pixels.
            # Project unit local east and north into tangent-plane components.
            u0=np.array([cos(d0),0,sin(d0)]);cx=np.array([0.,1,0]);cy=np.array([-sin(d0),0,cos(d0)])
            u=np.array([cos(d)*cos(da),cos(d)*sin(da),sin(d)])
            east=np.array([-sin(da),cos(da),0]);north=np.array([-sin(d)*cos(da),-sin(d)*sin(da),cos(d)])
            den=u@u0;j=np.array([[(v@t*den-(u@t)*(v@u0))/den**2 for v in [east,north]] for t in [cx,cy]])
            tan=radians(theta+90);ev,nv=np.linalg.solve(j,[-cos(tan),sin(tan)])
            effective=degrees(atan2(nv,-ev))
            stats=transform(stars,coords,ra,dec,effective)
            for label in ('axis','rotated_axis'):
                key=label+'_gaia_row'
                if key in stats:stats[label+'_source_id']=gaia[stats.pop(key)]['source_id']
            out.append(dict(object_id=r['object_id'],field=f,**stats))
    with (HERE/'gaia_alignment.csv').open('w',newline='') as dest:
        fields=['object_id','field','nearby_stars','coincident','aligned','rotated']+sorted({v for r in out for v in r if v.startswith(('axis_','rotated_axis_'))})
        w=csv.DictWriter(dest,fieldnames=fields);w.writeheader();w.writerows(out)
if __name__=='__main__':main()
