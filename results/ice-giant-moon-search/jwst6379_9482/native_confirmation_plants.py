"""Plant calibrated sources before the within-exposure RATEINTS median and apply compact-source vetting."""
import argparse,csv,json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from urllib.request import urlopen
import warnings
from astropy.io import fits
import numpy as np
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'jwst9482'));sys.path.insert(0,str(base/'jwst8975'))
from link_visits import geometry,search,trials
from native_plants import Detector
import positive_control as pc

FLAT_URL='https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_nircam_flat_0736.fits'


def make_frames(cube_path,native_path,meta,pixels,files,ee,flat):
    c=search.load_npz(cube_path)
    q=search.load_npz(native_path)
    rec=json.loads(Path(str(cube_path)+'.inputs.json').read_text())
    frames=[]
    for f,r in enumerate(rec):
        path=files/r['filename']
        d=Detector(path,r,c,ee)
        m=pc.metadata(meta,r['exposure_id'].removesuffix('_cal'))
        raw,x0,y0=pc.wide_images(pixels,m,subtract_template=False)
        h=fits.getheader(path,'SCI');p=fits.getheader(path)
        if p['R_FLAT'].split('/')[-1]!='jwst_nircam_flat_0736.fits':
            raise ValueError('unexpected flat calibration')
        jy_per_dn_s=h['PHOTMJSR']*h['PIXAR_SR']*1e6
        frames.append((d,raw,x0,y0,jy_per_dn_s,flat,q['frame']==f))
    return frames


def native_stack(frames,truth,claimed,k,flux):
    stamps=[]
    for d,raw,x0,y0,factor,flat,sel in frames:
        tx,ty=d.native(truth['x_path'][k,sel],truth['y_path'][k,sel])
        cx,cy=d.native(claimed['x'][k,sel],claimed['y'][k,sel])
        # Include the complete 49-pixel photometric region and broad-filter support.
        xmin=max(x0,int(np.floor(min(tx.min(),cx.min())-48)))
        xmax=min(x0+raw.shape[2],int(np.ceil(max(tx.max(),cx.max())+49)))
        ymin=max(y0,int(np.floor(min(ty.min(),cy.min())-48)))
        ymax=min(y0+raw.shape[1],int(np.ceil(max(ty.max(),cy.max())+49)))
        if xmax<=xmin or ymax<=ymin:
            stamps.extend([np.full((49,49),np.nan) for _ in cx]);continue
        vals=raw[:,ymin-y0:ymax-y0,xmin-x0:xmax-x0].copy()
        flat_cut=flat[ymin:ymax,xmin:xmax]
        for i in range(len(vals)):
            x8=np.rint((tx[i]-xmin)*8).astype(int)
            y8=np.rint((ty[i]-ymin)*8).astype(int)
            weights=np.ones(len(x8))/(len(x8)-1);weights[[0,-1]]*=.5
            for x,y,w in zip(x8,y8,weights,strict=True):
                patch=search.clipped_stamp(vals[i].shape,y//8-d.half,x//8-d.half,d.phase_grid[y%8,x%8])
                if patch is not None:
                    cut,psf=patch
                    vals[i][cut]+=psf*(w*flux/factor)*flat_cut[cut]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            vals-=np.nanmedian(vals,axis=0)
        stamps.extend(pc.patches_for_positions(vals,xmin,ymin,{'object':np.stack([cx,cy],axis=1)})['object'])
    arr=np.asarray(stamps);mid=len(arr)//2
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        co=[np.nanmedian(z,axis=0) for z in [arr,arr[:mid],arr[mid:]]]
    yy,xx=np.mgrid[-pc.RADIUS:pc.RADIUS+1,-pc.RADIUS:pc.RADIUS+1]
    ann=(xx**2+yy**2>=8**2)&(xx**2+yy**2<=20**2)
    ims=[]
    for im in co:
        m=np.nanmedian(im[ann]);s=1.4826*np.nanmedian(abs(im[ann]-m))
        ims.append((im-m)/s)
    inner=(xx**2+yy**2<=3**2)&np.isfinite(ims[0])
    if not inner.any():return [0.,0.,0.]
    ix=np.unravel_index(np.nanargmax(np.where(inner,ims[0],np.nan)),ims[0].shape)
    return [z[ix] for z in ims]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube1','cube2','map1','map2','meta1','pixels1','meta2','pixels2','files','ee','flat','replay','out']:p.add_argument(k,type=Path)
    a=p.parse_args()
    if not a.flat.exists():a.flat.write_bytes(urlopen(FLAT_URL,timeout=90).read())
    flat=fits.getdata(a.flat,'SCI')
    rows=list(csv.DictReader((a.replay/'plants.csv').open()))
    if not (a.replay/'report.json').exists():raise ValueError('complete discovery replay first')
    e0=json.loads((a.replay/'report.json').read_text())['epoch_jd_tdb']
    elem=np.zeros((len(rows),6));claimed=elem.copy()
    for k,r in enumerate(rows):
        elem[k,0],elem[k,-1]=float(r['a_km']),float(r['phase_epoch1'])
        # Unrecovered trials remain failures; dummy coordinates are never measured.
        claimed[k]=elem[k]
        if r['recovered']=='True':
            claimed[k,0],claimed[k,-1]=float(r['recovered_a_km']),float(r['recovered_phase_deg'])
    a.out.mkdir(exist_ok=False,parents=True)
    results=[[] for _ in rows]
    with TemporaryDirectory() as td:
        tmp=Path(td)
        for path,mp,meta,pixels in [(a.cube1,a.map1,a.meta1,a.pixels1),(a.cube2,a.map2,a.meta2,a.pixels2)]:
            truth=trials(tmp/'e.csv',mp,tmp/'q.npz',elem,e0)
            fitted=trials(tmp/'e.csv',mp,tmp/'q.npz',claimed,e0)
            frames=make_frames(path,mp,meta,pixels,a.files,a.ee,flat)
            for k,r in enumerate(rows):
                stats=native_stack(frames,truth,fitted,k,float(r['flux'])) if r['recovered']=='True' else [0.,0.,0.]
                results[k].extend(stats)
    cols=['flux','plant_index','cal_recovered','recovered','status',
          'feb_peak','feb_first','feb_second','oct_peak','oct_first','oct_second']
    final=[]
    with (a.out/'native_dispositions.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        for r,v in zip(rows,results,strict=True):
            z=np.asarray(v).reshape(2,3)
            passed=bool((z[:,0]>=5).all() and (z[:,1:]>=2).all())
            row=dict(zip(cols,[float(r['flux']),int(r['plant_index']),r['recovered']=='True',passed,'recovered' if passed else 'missed',*v],strict=True))
            w.writerow(row);final.append(row)
    summary=dict(native_flat_url=FLAT_URL,flat_sha256=search.file_hash(a.flat),code_sha256=search.file_hash(Path(__file__)),
                 local_contrast=5,each_chronological_half_at_same_peak=2,local_peak_search_radius_px=3,
                 calibrated_source_noise='Fixed observed noise and DQ with model flux added to native RATEINTS.',
                 joint_recovery=search.summarize_plants(final),replay_sha256=search.file_hash(a.replay/'plants.csv'))
    (a.out/'report.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
