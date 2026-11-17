"""Independent interpolation/count replay of the same two authorized PUNCH images."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import scipy
from scipy.ndimage import map_coordinates
import astropy
from astropy.io import fits
from astropy.wcs import WCS

base=Path('D:/Open-Space/campaigns/punch-root-independent-16f73d35-20261117')
base.mkdir(exist_ok=False)
campaign=Path('D:/Open-Space/campaigns/punch-development-16f73d35-20261117')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(campaign/'manifest.json')=='4c26bcca66260877e802078abd75281fcf8ce6cd4b9778b3983de4697033b91e'
assert sha(campaign/'authorization.json')=='dd3190f149e8c30a06f47c78b16bccc1e39c61eb9edfda901456a7f61a932d66'
primary_sha=sha(campaign/'development.json')
assert primary_sha=='ae5b77aeb30e898a72a6f160b9fcc5ad7a05682b9ae2a5676ba1982d165930e6'
receipt=json.loads((campaign/'execution.json').read_bytes())
assert receipt['output_sha256']==primary_sha and receipt['exit_code']==0
assert sha(campaign/'run.log')==receipt['log_sha256']
assert sha(campaign/'command.json')==receipt['command_sha256']
m=json.loads((campaign/'manifest.json').read_bytes())
for row in m['inputs']:
    assert sha(Path(row['path']))==row['sha256']
preflight={'scope':'Same two first-development mosaics only; no prospective input',
           'code_sha256':sha(Path(__file__)),'primary_output_sha256':primary_sha,
           'inputs':m['inputs'],'numpy':np.__version__,'astropy':astropy.__version__,
           'scipy':scipy.__version__,'python':sys.version}
(base/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
images=[];masks=[];world=[];coverage=[]
for row in m['inputs']:
    with fits.open(row['path'],memmap=False) as h:
        image=np.array(h['PRIMARY DATA ARRAY'].data,dtype=np.float64)
        packed=np.array(h['UNCERTAINTY ARRAY'].data,dtype=np.float64)
        world.append(WCS(h['PRIMARY DATA ARRAY'].header,naxis=2))
    with np.errstate(divide='ignore',invalid='ignore'):
        sigma=np.divide(image,packed)
    valid=np.isfinite(image)&np.isfinite(packed)&(packed!=0)&np.isfinite(sigma)&(sigma>0)
    coverage.append({'total':int(image.size),'finite_image':int(np.isfinite(image).sum()),
                     'finite_packed':int(np.isfinite(packed).sum()),'zero_packed':int((packed==0).sum()),
                     'finite_positive_decoded_sigma':int((np.isfinite(sigma)&(sigma>0)).sum()),
                     'valid':int(valid.sum())})
    images.append(image);masks.append(valid)
a,b=images;ma,mb=masks
height,width=b.shape
registered=np.full(a.shape,np.nan)
for start in range(0,a.shape[0],96):
    y,x=np.mgrid[start:min(start+96,a.shape[0]),:a.shape[1]]
    longitude,latitude=world[0].all_pix2world(x,y,0)
    qx,qy=world[1].all_world2pix(longitude,latitude,0)
    inside=np.isfinite(qx)&np.isfinite(qy)&(qx>=0)&(qy>=0)&(qx<width-1)&(qy<height-1)
    ix=np.floor(np.where(inside,qx,0)).astype(int)
    iy=np.floor(np.where(inside,qy,0)).astype(int)
    accepted=inside&mb[iy,ix]&mb[iy+1,ix]&mb[iy,ix+1]&mb[iy+1,ix+1]
    coords=np.stack((np.where(inside,qy,0),np.where(inside,qx,0)))
    interpolation=map_coordinates(b,coords,order=1,mode='constant',cval=np.nan,prefilter=False)
    registered[start:start+y.shape[0]]=np.where(accepted,interpolation,np.nan)
valid=ma&np.isfinite(registered)
delta=registered-a
values=delta[valid]
summary=json.loads((campaign/'development.json').read_bytes())
assert summary['status']=='completed'
expected=summary['result']
assert [r['valid'] for r in coverage]==expected['native_valid_pixels']
assert int(valid.sum())==expected['eligible_pixels']
stats=[float(np.median(values)),*map(float,np.percentile(values,[5,95]))]
target=[expected['difference_median'],*expected['difference_p05_p95']]
np.testing.assert_allclose(stats,target,rtol=1e-9,atol=1e-25)
grid=[];errors=[];finite_blocks=0
for i,ys in enumerate(np.array_split(np.arange(a.shape[0]),64)):
    row=[]
    for j,xs in enumerate(np.array_split(np.arange(a.shape[1]),64)):
        block=np.ix_(ys,xs);v=delta[block][valid[block]]
        found=float(np.mean(v)) if v.size else None
        recorded=expected['difference_grid_64'][i][j]
        assert (found is None)==(recorded is None)
        if found is not None:
            errors.append(abs(found-recorded));finite_blocks+=1
            assert abs(found-recorded)<=1e-25+abs(recorded)*1e-9
        row.append(found)
    grid.append(row)
report={'status':'passed','same_observations':True,'coverage':coverage,'common_eligible':int(valid.sum()),
        'descriptive_summary':stats,'finite_grid_blocks':finite_blocks,'max_grid_numeric_delta':max(errors),
        'grid':grid,'primary_output_sha256':primary_sha,
        'limits':'Independent interpolation implementation on the same two mosaics, not independent observations. No p-value, moving front, calibrated error model, validated quality mask or physical interpretation; later prospective mosaic remains unopened.'}
(base/'independent_check.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['status','coverage','common_eligible','finite_grid_blocks','max_grid_numeric_delta']}))
