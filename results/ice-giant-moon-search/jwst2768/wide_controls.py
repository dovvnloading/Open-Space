"""Check the WCS registration of the four outer bright moons using their native visit mosaics."""
import argparse,csv
from pathlib import Path
import sys
import numpy as np
from scipy import ndimage
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from geometry import sky_xy
from ephemeris import interpolate
from search import load_npz


def run(files,eph,out):
    rows=[]
    for file in sorted(files):
        m=load_npz(file);_,_,ids,obs=interpolate(eph,m['mjd'])
        for body in [701,702,703,704,705]:
            k=list(ids).index(body)
            off=sky_xy(obs[k,:,0],obs[k,:,1],obs[0,:,0],obs[0,:,1],1.,(0.,0.))
            off[:,0]*=-1
            path=off@m['M']+m['center']; x,y=path.mean(0)
            yy,xx=np.indices((33,33));left=int(round(x))-16;top=int(round(y))-16
            if left<0 or top<0 or left+33>m['z'].shape[1] or top+33>m['z'].shape[0]:continue
            sl=np.s_[top:top+33,left:left+33];good=m['raw_good'][sl];z=m['z'][sl]
            rho=np.hypot(xx+left-x,yy+top-y)
            ring=good&(rho>4)&(rho<12); high=(ring&(z>=5)).sum()
            # Use capped positive weights around the predicted sky position.
            # A saturated center may have no valid samples, so this is a
            # broad-halo registration check, not precise point astrometry.
            w=np.where(good&(rho<10),np.clip(z,0,30),0)
            if w.sum():
                mx=(w*(xx+left)).sum()/w.sum();my=(w*(yy+top)).sum()/w.sum()
            else:mx=my=np.nan
            rows.append([str(m['visit']),body,x,y,int(ring.sum()),int(high),np.nan if not ring.any() else np.percentile(z[ring],90),
                         mx,my,np.hypot(mx-x,my-y),m['pixel_scale']])
    with open(out,'w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['visit','body','predicted_x','predicted_y','usable_ring_pixels','ring_pixels_score_ge_5',
                                     'ring_p90_score','halo_moment_x','halo_moment_y','halo_centroid_separation_px','scale_arcsec']);wr.writerows(rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('ephemeris');p.add_argument('output');p.add_argument('maps',nargs='+')
    a=p.parse_args();run(a.maps,a.ephemeris,a.output)
