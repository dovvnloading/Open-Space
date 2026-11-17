"""Fixed finite surface footprints; header geometry only, no image access."""
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from sunpy.coordinates import HeliographicStonyhurst,Helioprojective,propagate_with_solar_surface
from geometry import to_surface,observer,mu,dust_center_zero_based

GRID=(-.4,-.2,0.,.2,.4)

def patch_lattice():
    # Exact finite footprint: 21x21 tangent-plane nodes inside a 0.02-radian
    # solar-surface circle, including its boundary. All interpolation stencil
    # pixels are separately checked, so no untested node contributes.
    z=np.linspace(-.02,.02,21);x,y=np.meshgrid(z,z);keep=x*x+y*y<=.02**2+1e-15
    return x[keep],y[keep]

def detector_fixed(reference, target, reference_mapping):
    """Same physical detector-coordinate nodes, not the tracked solar patch."""
    x=np.asarray(reference_mapping['x'])+reference['PXBEG1']-target['PXBEG1']
    y=np.asarray(reference_mapping['y'])+reference['PXBEG2']-target['PXBEG2']
    c=to_surface(target,x,y);m=np.asarray(mu(target,c),float)
    finite=np.isfinite(x)&np.isfinite(y)&np.isfinite(m)
    inside=finite&(x>=0)&(y>=0)&(x<target['NAXIS1']-1)&(y<target['NAXIS2']-1)
    cx,cy=dust_center_zero_based(target);nodust=(np.abs(x-cx)>21)|(np.abs(y-cy)>21)
    valid=inside&(m>=.65)&nodust
    return {'x':x,'y':y,'mu':m,'eligible':bool(np.all(valid)),'invalid_nodes':int((~valid).sum())}


def footprints(reference,headers):
    theta=np.rad2deg(np.arcsin(reference['RSUN_REF']/reference['DSUN_OBS']))
    w=WCS(reference,naxis=2);dx,dy=patch_lattice();result=[];centers=[]
    for i,x in enumerate(GRID):
        for j,y in enumerate(GRID):
            px,py=w.all_world2pix(x*theta,y*theta,0)
            c=to_surface(reference,px,py)
            base=SkyCoord(c.transform_to(HeliographicStonyhurst(obstime=reference['DATE-AVG'])))
            if any(base.separation(previous).to_value(u.rad)<=.04 for previous in centers):raise ValueError('surface patches overlap')
            centers.append(base)
            local=base.skyoffset_frame()
            points=SkyCoord(lon=dx*u.rad,lat=dy*u.rad,distance=reference['RSUN_REF']*u.m,frame=local).transform_to(base.frame)
            maps=[];why=[]
            for h in headers:
                target=Helioprojective(observer=observer(h),obstime=h['DATE-AVG'],rsun=h['RSUN_REF']*u.m)
                with propagate_with_solar_surface(rotation_model='howard'):v=points.transform_to(target)
                visible=np.asarray(v.is_visible(),bool);m=np.asarray(mu(h,v),float)
                xp,yp=WCS(h,naxis=2).all_world2pix(v.Tx.to_value(u.deg),v.Ty.to_value(u.deg),0)
                xp=np.asarray(xp,float);yp=np.asarray(yp,float)
                finite=np.isfinite(xp)&np.isfinite(yp)&np.isfinite(m)
                inside=finite&(xp>=0)&(yp>=0)&(xp<h['NAXIS1']-1)&(yp<h['NAXIS2']-1)
                # Conservative 40x40 square around the release team's moving
                # one-based 30x30 dust center, plus bilinear-stencil margin.
                cx,cy=dust_center_zero_based(h)
                nodust=(np.abs(xp-cx)>21)|(np.abs(yp-cy)>21)
                valid=visible&inside&(m>=.65)&nodust
                if not np.all(valid):why.append({'filename':h['FILENAME'],'visible_fail':int((~visible).sum()),'mu_fail':int((~(m>=.65)).sum()),'edge_fail':int((~inside).sum()),'dust_fail':int((~nodust).sum())})
                maps.append({'x':xp,'y':yp,'mu':m})
            result.append({'patch_id':f'g{i}{j}','grid_fraction':[x,y],'eligible':not why,'reasons':why,'maps':maps})
    return result
