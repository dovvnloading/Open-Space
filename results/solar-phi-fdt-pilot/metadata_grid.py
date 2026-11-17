"""Header-only grid-center feasibility; footprint and pixels remain untested."""
import argparse
import json
import numpy as np
from astropy.wcs import WCS
from geometry import read_header,to_surface,mu,transport

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('development');p.add_argument('repeat_a');p.add_argument('repeat_b');a=p.parse_args()
    h=[read_header(x) for x in (a.development,a.repeat_a,a.repeat_b)]
    theta=np.rad2deg(np.arcsin(h[0]['RSUN_REF']/h[0]['DSUN_OBS']))
    w=WCS(h[0],naxis=2);rows=[]
    for x in (-.4,-.2,0.,.2,.4):
        for y in (-.4,-.2,0.,.2,.4):
            px,py=w.all_world2pix(x*theta,y*theta,0)
            c=to_surface(h[0],px,py); mus=[float(mu(h[0],c))];visible=[]
            for target in h[1:]:
                xp,yp,v=transport(h[0],target,px,py);visible.append(bool(v))
                if v:
                    cp=to_surface(target,xp,yp);mus.append(float(mu(target,cp)))
                else:mus.append(None)
            rows.append({'x_fraction':x,'y_fraction':y,'mu_centers':mus,'visible':visible,'passes_center_only':all(visible) and all(v>=.65 for v in mus if v is not None)})
    print(json.dumps({'reference':h[0]['FILENAME'],'center_only':True,'eligible_centers':sum(z['passes_center_only'] for z in rows),'grid':rows},indent=2))
