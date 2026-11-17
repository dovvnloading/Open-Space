"""Evaluate registered tangent-plane residuals and orbit-matched control positions."""
import argparse
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from scipy import ndimage as ndi
from link_visits import geometry,search,make_cube,trials,evaluate
ROOT=Path(__file__).resolve().parent


def save(path, records):
    with path.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(records[0]));w.writeheader();w.writerows(records)


def main(args):
    best=max(csv.DictReader((ROOT/'grid/refined.csv').open()),key=lambda r:min(float(r['deep_score']),float(r['confirmation_score'])))
    a,phi=float(best['a_km']),float(best['phase_deg'])
    frames=[]; comparisons=[];summary={}
    with TemporaryDirectory() as tmp:
        d=Path(tmp)
        geometry.generate(args.deep,d/'epoch.npz',grid=[a,a,1,360])
        epoch=float(np.load(d/'epoch.npz')['epoch_jd_tdb'])
        for label,path,deep in [('F150W2',args.deep,True),('F162M',args.confirm,False)]:
            c=make_cube(path,args.ee,deep)
            psfs=search.validate_cube(c,min_visits=1);cuts=search.Cuts(psf_support=.8);maps=search.prepare(c,psfs,cuts)
            offsets=np.r_[0.,np.arange(-32,0,1.),np.arange(1,33,1.)]
            elem=np.zeros((len(offsets),6));elem[:,0]=a;elem[:,-1]=(phi+offsets)%360
            bank=trials(d/'t.csv',path,d/'b.npz',elem,epoch)
            n,i=evaluate(c,bank,psfs,maps,cuts)
            z=search.snr(n,i);rawflux=np.divide(n,i,out=np.zeros_like(n),where=i>0)
            # Compare tracks radially and perpendicular to the exposure trail.
            # Local medians use only placements on the same equatorial orbit.
            far=(abs(offsets)>=5)&(abs(offsets)<=20)
            centered=n[0].copy(); pm=np.stack([bank['x_path'][0],bank['y_path'][0]],axis=-1)
            for j,fid in enumerate(c['exposure_id']):
                valid=far&(i[:,j]>0)
                sky=np.median(rawflux[valid,j]) if valid.any() else 0.
                centered[j]-=sky*i[0,j]
                # Minimum over all sample pairs; the entire exposure is integrated.
                distances=np.linalg.norm(c['known_path'][:,j,:,None,:]-pm[None,j,None,:,:],axis=-1)
                dmin=distances.min(axis=(1,2));idx=dmin.argmin()
                frames.append(dict(band=label,frame=str(fid),visit=str(c['visit'][j]),supported=int(i[0,j]>0),
                            score=round(z[0,j],3),local_flux_uJy=round(sky*1e6,5),target_flux_uJy=round(rawflux[0,j]*1e6,5),
                            nearest_trail_body=int(c['known_id'][idx]),min_all_pairs_px=round(dmin[idx],3),
                            same_radius_controls=int(valid.sum())))
            for v in np.unique(c['visit']):
                use=c['visit']==v
                raw=search.snr(n[:,use].sum(1),i[:,use].sum(1))
                cnt=(i[:,use]>0).sum(1);equal=(cnt==cnt[0])&(offsets!=0)
                # Matched translations for the same trail times and orientation.
                trans=np.array([(dx,dy) for dx in range(-4,5) for dy in range(-4,5)])
                shifted={key:np.repeat(value[:1],len(trans),axis=0) for key,value in bank.items()
                         if isinstance(value,np.ndarray) and value.ndim and value.shape[0]==len(offsets)}
                for key in ['x','x_path']: shifted[key]+=trans[:,0].reshape((-1,)+(1,)*(shifted[key].ndim-1))
                for key in ['y','y_path']: shifted[key]+=trans[:,1].reshape((-1,)+(1,)*(shifted[key].ndim-1))
                sn,si=evaluate(c,shifted,psfs,maps,cuts)
                zs=search.snr(sn[:,use].sum(1),si[:,use].sum(1));same=(si[:,use]>0).sum(1)==cnt[0]
                peak=np.argmax(np.where(same,zs,-1e10))
                comparisons.append(dict(band=label,visit=str(v),used=int(cnt[0]),rank_score=round(float(raw[0]),3),
                    flank_corrected_score=round(float(search.snr(centered[use].sum(keepdims=True),i[0,use].sum(keepdims=True))[0]),3),
                    phase_controls=int(equal.sum()),phase_controls_ge=int((equal&(raw>=raw[0])).sum()),
                    translated_peak_dx=int(trans[peak,0]),translated_peak_dy=int(trans[peak,1]),translated_score=round(float(zs[peak]),3)))
            summary[label]=dict(cube_sha256=search.file_hash(path),radius_km=a,phase_deg=phi,emission_epoch_tdb_jd=epoch)
    save(ROOT/'independent_frames.csv',frames);save(ROOT/'independent_visits.csv',comparisons)
    (ROOT/'independent_vet.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('deep','confirm','ee'):ap.add_argument(k,type=Path)
    main(ap.parse_args())
