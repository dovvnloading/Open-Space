"""Independent native CF unpacking checks of six authorized G19 granules only."""
from pathlib import Path
import collections,datetime,hashlib,json,sys
import h5py,numpy as np
root=Path('D:/Open-Space/campaigns/goes19-characterization-4d8225ab-20261117')
out=Path('D:/Open-Space/campaigns/goes19-root-independent-4d8225ab-20261117')
out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(root/'manifest.json')=='df2d900ab9fb7ff348e985207cf15a20ab34745a9681fc7bb5fc980cf16f299f'
assert sha(root/'authorization.json')=='92ec544c05d67209bc920152dd081f7bc04f9ff05334874e97b4a5c3c6dc9e09'
assert sha(root/'characterization.json')=='500c519ad637ae07d0916ab964056a3da2abb18c3599e6e6bebf9c34263cf268'
manifest=json.loads((root/'manifest.json').read_bytes())
primary=json.loads((root/'characterization.json').read_bytes())
assert primary['status']=='completed' and len(primary['files'])==len(manifest['inputs'])==6
for row in manifest['inputs']:
    assert row['satellite']==19 and sha(row['path'])==row['sha256']
preflight={'code_sha256':sha(__file__),'inputs':manifest['inputs'],'primary_sha256':sha(root/'characterization.json'),
           'python':sys.version,'numpy':np.__version__,'h5py':h5py.__version__,
           'scope':'Same six authorized G19 files only; no G18 event arrays.'}
(out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')

def decode(dataset):
    raw=dataset[()]
    good=np.ones(raw.shape,dtype=bool)
    if '_FillValue' in dataset.attrs:
        fill=np.asarray(dataset.attrs['_FillValue']).reshape(-1)[0]
        good &= raw!=fill
    unsigned=str(dataset.attrs.get('_Unsigned','')).lower() in ('true',"b'true'")
    numeric=raw.astype(np.int64) if raw.dtype.kind in 'iu' else raw.astype(np.float64)
    if unsigned:
        numeric=np.bitwise_and(numeric,(1<<(8*raw.dtype.itemsize))-1)
    return numeric,good

results=[];maximum_time_delta=0.0
for index,(row,record) in enumerate(zip(manifest['inputs'],primary['files'])):
    assert row['sha256']==record['sha256']
    with h5py.File(row['path'],'r') as h:
        time,good=decode(h['group_time_offset'])
        scale=float(np.asarray(h['group_time_offset'].attrs['scale_factor']).reshape(-1)[0])
        offset=float(np.asarray(h['group_time_offset'].attrs['add_offset']).reshape(-1)[0])
        elapsed=time*scale+offset+index*20
        identifiers,good_id=decode(h['group_id'])
        parent,good_parent=decode(h['group_parent_flash_id'])
        quality,good_quality=decode(h['group_quality_flag'])
        lat=h['group_lat'][()].astype(np.float64);lon=h['group_lon'][()].astype(np.float64)
    assert good.all() and good_id.all() and good_parent.all() and good_quality.all()
    assert len(record['rows'])==len(elapsed)
    assert [r['row_index'] for r in record['rows']]==list(range(len(elapsed)))
    assert all(r['file_sha256']==row['sha256'] for r in record['rows'])
    numeric_time=np.array([r['elapsed_utc_s'] for r in record['rows']])
    delta=float(np.max(np.abs(numeric_time-elapsed)))
    maximum_time_delta=max(maximum_time_delta,delta)
    np.testing.assert_allclose(numeric_time,elapsed,rtol=0,atol=1e-12)
    for field,native in [('group_id',identifiers),('parent_flash_id',parent),('quality',quality),('lat_deg',lat),('lon_deg',lon)]:
        np.testing.assert_array_equal(np.array([r[field] for r in record['rows']]),native)
    inside=(elapsed>=0)&(elapsed<120)
    np.testing.assert_array_equal(np.array([r['inside_interval'] for r in record['rows']]),inside)
    counts={'attempted':len(elapsed),'inside':int(inside.sum()),'outside':int((~inside).sum()),
            'good_qc':int((quality==0).sum()),'degraded_qc':int(np.isin(quality,[1,3,5]).sum())}
    for key,value in counts.items():assert record['counts'][key]==value
    assert np.isfinite(lat).all() and np.isfinite(lon).all() and (np.abs(lat)<=90).all() and (np.abs(lon)<=180).all()
    assert all(record['counts'][key]==0 for key in ['invalid_time','unknown_or_fill_qc','invalid_positions','invalid_id_or_parent']) if 'invalid_positions' in record['counts'] else all(record['counts'][key]==0 for key in ['invalid_time','unknown_or_fill_qc','invalid_position','invalid_id_or_parent'])
    results.append({'file_sha256':row['sha256'],'counts':counts,'time_numeric_delta':delta})
answer={'status':'native_rows_match','files':results,'totals':{k:sum(r['counts'][k] for r in results) for k in results[0]['counts']},
        'max_time_numeric_delta':maximum_time_delta,'limits':'Same six files and file-local rows, not independent physical flashes. No deduplication, candidate selection, height, p-value, duty-cycle completeness or nondetection inference. G18 event arrays remain sealed.'}
(out/'independent_check.json').write_text(json.dumps(answer,indent=2)+'\n')
print(json.dumps({'status':answer['status'],'totals':answer['totals'],'max_time_numeric_delta':maximum_time_delta,'report':str(out/'independent_check.json')}))
