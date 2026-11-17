"""Force 2008 NICI/blue H-band apertures on both astrometric hypotheses in each exposure.
Diagnostic on detector images, with matched same-radius halo controls; the independently
published full-ADI 95% limits and their PA/coverage caveats are in NICI_HD38206_loci.csv.
"""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi
from scipy.optimize import least_squares
def main(parent,cache,out):
    sys.path.insert(0,str(parent.resolve()/'vetting'))
    from archival_images import test_locus
    inp=list(csv.DictReader((out/'gemini_inputs.csv').open()));loci=list(csv.DictReader((out/'NICI_HD38206_loci.csv').open()))
    results=[];meta=[]
    for row in inp:
        path=cache/row['file']
        with path.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==row['compressed_SHA256']
        with fits.open(path) as f:
            p=f[0].header;assert p['FILTER_B'].startswith('H_')
            hd=next(h for h in f[1:] if h.header['CHANNEL']=='BLUE')
            h=hd.header;image=hd.data.astype(float);x,y=np.rint([h['CRPIX1']-1,h['CRPIX2']-1]).astype(int)
            # Fit the centre of point symmetry of the core/halo of the occulted host.
            # Background shifts cancel by comparing each sample with its rotated counterpart.
            yy,xx=np.mgrid[y-48:y+49:2,x-48:x+49:2]
            smooth=ndi.median_filter(image,size=3)
            sky=np.nanmedian(smooth);sub=smooth[yy,xx];weight=np.maximum(sub-sky,0)
            init=np.array([np.sum(xx*weight),np.sum(yy*weight)])/np.sum(weight)
            def diff(v):
                mirror=ndi.map_coordinates(smooth,[2*v[1]-yy,2*v[0]-xx],order=1,mode='nearest')
                return (sub-mirror).ravel()
            center=least_squares(diff,init,loss='soft_l1',f_scale=200).x
            mat=3600*np.array([[h['CD1_1'],h['CD1_2']],[h['CD2_1'],h['CD2_2']]])
            scl=abs(np.linalg.det(mat))**.5;assert .016<scl<.020
            meta.append(dict(input=row['file'],host_x=center[0],host_y=center[1],scale=scl,rot=np.rad2deg(np.arctan2(*mat[:,1]))))
            # Slow-varying halo and static detector bias locally cancel in test_locus.
            for l in loci:
                delta=np.array([float(l['east_arcsec']),float(l['north_arcsec'])])
                pred=center+np.linalg.solve(mat,delta)
                r=test_locus(image,None,pred,center,.075,scl)
                results.append(dict(input=row['file'],visit=l['visit'],x=l['x'],y=l['y'],hypothesis=l['hypothesis'],pixel_x=pred[0],pixel_y=pred[1],
                                    coverage=r['coverage'],raw_halo_score=r.get('score',''),n_control=r.get('n_control','')))
    def write(name,rows):
        with (out/name).open('w') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
    write('gemini_image_checks.csv',results);write('gemini_centers.csv',meta)
    summary=[]
    for l in loci:
        sel=[r for r in results if all(r[k]==l[k] for k in ['x','y','visit','hypothesis']) and r['coverage']=='measured']
        scores=[r['raw_halo_score'] for r in sel]
        summary.append(dict(visit=l['visit'],x=l['x'],y=l['y'],hypothesis=l['hypothesis'],measured=len(sel),median_raw_score=np.median(scores),min_score=min(scores),max_score=max(scores),N_gt_5=sum(v>=5 for v in scores)))
    write('gemini_locus_summary.csv',summary)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('cache',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.parent,a.cache,a.out)
