"""Source adapter, authorization-gated. No validation-sensor reader exists."""
import hashlib,json,re,datetime as dt
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def validate(manifest,auth,manifest_digest,auth_digest,paths):
    if sha(manifest)!=manifest_digest or sha(auth)!=auth_digest: raise ValueError('manifest/auth digest')
    m=json.loads(Path(manifest).read_text()); a=json.loads(Path(auth).read_text())
    if m.get('phase')!='source-development' or a.get('phase')!='source-development' or a.get('manifest_sha256')!=manifest_digest: raise ValueError('wrong phase')
    if m.get('reader_sha256')!=sha(__file__): raise ValueError('reader digest')
    if len(paths)!=6 or len(m.get('inputs',[]))!=6: raise ValueError('exactly six source inputs')
    for index,(path,row) in enumerate(zip(paths,m['inputs'])):
        p=Path(path)
        if not p.is_absolute() or p.drive.upper()!='D:' or str(p.resolve())!=row['path'] or row['satellite']!=19 or '_G19_' not in p.name or sha(p)!=row['sha256']: raise ValueError('wrong source input')
        second=index*20; stamp=f'{second//60:02d}{second%60:02d}0'
        if '_s202632000'+stamp+'_e' not in p.name: raise ValueError('unexpected time sequence')
    return m

def decode_groups(handle):
    import numpy as np
    required=('group_time_offset','group_lat','group_lon','group_quality_flag','group_energy','group_parent_flash_id','group_id')
    if any(k not in handle for k in required): raise ValueError('missing group field')
    d=handle['group_time_offset']; attrs=d.attrs
    if attrs.get('_Unsigned',b'') not in (b'true','true') or d.dtype.kind!='i' or d.dtype.itemsize!=2: raise ValueError('time packing schema')
    units=attrs.get('units',b''); units=units.decode() if isinstance(units,bytes) else units
    if not units.startswith('seconds since '): raise ValueError('time units')
    origin=dt.datetime.fromisoformat(units.removeprefix('seconds since ')).replace(tzinfo=dt.timezone.utc)
    base=dt.datetime(2026,11,16,tzinfo=dt.timezone.utc)
    if not (0 <= (origin-base).total_seconds() < 120): raise ValueError('time origin outside cohort')
    scale=float(np.asarray(attrs['scale_factor']).item()); offset=float(np.asarray(attrs['add_offset']).item())
    if not np.isfinite(scale) or scale<=0 or not np.isfinite(offset): raise ValueError('invalid time packing')
    packed=np.asarray(d[:]); t=packed.view(np.uint16).astype(float)*scale+offset
    lat=np.asarray(handle['group_lat'][:],float);lon=np.asarray(handle['group_lon'][:],float); quality=np.asarray(handle['group_quality_flag'][:]);
    if any(len(z)!=len(t) for z in (lat,lon,quality)) or not np.isin(quality,[0,1,3,5]).all(): raise ValueError('group shape/quality')
    # Data retain 20-second boundary margin for attribution; cohort crop is applied later.
    elapsed=t+(origin-base).total_seconds(); ids=np.asarray(handle['group_id'][:]); parent=np.asarray(handle['group_parent_flash_id'][:])
    if len(ids)!=len(t) or len(parent)!=len(t): raise ValueError('identifier shape')
    keep=np.isfinite(elapsed)&(elapsed>=0)&(elapsed<120)
    return {'elapsed':elapsed[keep],'lat':lat[keep],'lon':lon[keep],'quality':quality[keep],'group_id':ids[keep],'parent':parent[keep],'attempted':len(t),'inside_interval':int(keep.sum())}

def read_source(manifest,auth,manifest_digest,auth_digest,paths):
    validate(manifest,auth,manifest_digest,auth_digest,paths)
    import h5py
    result=[]
    for path in paths:
        with h5py.File(path,'r') as f: result.append(decode_groups(f))
    return result
