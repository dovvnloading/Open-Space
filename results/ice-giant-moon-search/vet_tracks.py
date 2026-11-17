"""Use equal-radius flanking trails and coordinate checks on screened-in templates."""
import argparse
import csv
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from geometry import ARCSEC, generate, sky_xy
from prepare import inverse_sky
from search import Cuts, load_npz, prepare, sample_paths, score, validate_cube


def run(cube, data_dir, centers, candidates, output):
    import tempfile
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    c=load_npz(cube);cut=Cuts()
    m=prepare(c,validate_cube(c),cut)
    with open(centers) as f:
        centers={r['exposure_id']:(float(r['x']),float(r['y'])) for r in csv.DictReader(f)}
    with open(candidates) as f:
        groups={}
        for r in csv.DictReader(f):
            groups.setdefault(float(r['a_km']),[]).append(float(r['phase_deg']))
    offset=np.array([-15,-12,-10,-8,-6,6,8,10,12,15])
    summary=[]; frames=[]
    with tempfile.TemporaryDirectory() as dst:
        for a,phases in groups.items():
            # One-degree grid matches the parent screen and all prescribed flank samples.
            generate(cube,Path(dst)/'t.npz',grid=(a,a,1,1))
            t=load_npz(Path(dst)/'t.npz')
            if any(p!=int(p) for p in phases):
                raise ValueError("use the originating one-degree phase grid")
            n,w=sample_paths(m,t['x_path'],t['y_path'],cut)
            flux=np.divide(n,w,out=np.full_like(n,np.nan),where=w>0)
            for phase in phases:
                i=int(phase);near=(i+offset)%360
                local=np.nanmedian(flux[near],axis=0)
                ok=np.isfinite(local)&(np.isfinite(flux[near]).sum(0)>=2)
                wn=np.where(ok,w[i],0)
                nn=np.where(ok,n[i]-local*w[i],0)
                fit=score(nn[None,:],wn[None,:],c['visit'],cut)
                x,y=t['x'][i],t['y'][i]
                e,north=-(x-c['center_xy'][0])*c['scale_arcsec']/ARCSEC,(y-c['center_xy'][1])*c['scale_arcsec']/ARCSEC
                k=c['path_mjd'].shape[1]//2
                ra,dec=[],[]
                for j,exp in enumerate(c['exposure_id']):
                    r,d=inverse_sky(c['path_ra'][j,k],c['path_dec'][j,k],e[j],north[j]);ra.append(r);dec.append(d)
                    hdr=fits.getheader(Path(data_dir)/(exp+'.fits'),'SCI')
                    wcs=WCS(hdr,naxis=2)
                    cra,cdec=wcs.all_pix2world(*centers[exp],0)
                    nr,nd=inverse_sky(cra,cdec,e[j],north[j]);dx,dy=wcs.all_world2pix(nr,nd,0)
                    pos=c['known_path'][:,j]
                    dist=np.hypot(pos[:,: ,0]-t['x_path'][i,j],pos[:,:,1]-t['y_path'][i,j]).min(1)
                    nearest=dist.argmin()
                    frames.append([a,phase,exp,int(c['known_id'][nearest]),float(dist[nearest]),float(dx),float(dy),
                                   float(r),float(d),int(np.isfinite(flux[near,j]).sum()),flux[i,j]*1e6,local[j]*1e6])
                sid=sky_xy(np.array(ra),np.array(dec),ra[0],dec[0],1.,(0.,0.))
                pair=np.linalg.norm(sid[:,None]-sid[None,:],axis=2)
                summary.append([a,phase,float(pair.max()),*fit['visit_snr'][0],fit['snr'][0],
                                *fit['n'][0],bool(fit['rank'][0]>=1)])
    for name,header,rows in [
        ('flanking_trails.csv',['a_km','phase_deg','sidereal_span_arcsec','visit1_control_score','visit2_control_score',
                                'pooled_control_score','n_visit1','n_visit2','pass_after_flanking_subtraction'],summary),
        ('frame_checks.csv',['a_km','phase_deg','exposure_id','nearest_known_id','min_trail_distance_px','native_x','native_y',
                             'mid_ra_deg','mid_dec_deg','n_usable_flanks','raw_fit_uJy','median_flank_uJy'],frames)]:
        with (out/name).open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(header);writer.writerows(rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['cube','data_dir','centers','candidates','output']:
        p.add_argument(arg)
    a=p.parse_args()
    run(a.cube,a.data_dir,a.centers,a.candidates,a.output)
