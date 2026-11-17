"""Six-file G19 row characterization; does not expose a GOES-18 reader."""
import argparse,datetime as dt,hashlib,json,math,sys
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def path(p):
 p=Path(p)
 if not p.is_absolute() or p.drive.upper()!='D:': raise ValueError('absolute D: path required')
 return str(p.resolve())
def authorize(manifest,auth,manifest_hash,auth_hash,paths,output):
 if sha(manifest)!=manifest_hash or sha(auth)!=auth_hash: raise ValueError('digest mismatch')
 m=json.loads(Path(manifest).read_text());a=json.loads(Path(auth).read_text()); out=path(output)
 if m.get('phase')!='g19-characterization' or a.get('phase')!='g19-characterization' or a.get('manifest_sha256')!=manifest_hash or a.get('output_path')!=out or m.get('output_path')!=out: raise ValueError('unauthorized phase/output')
 if Path(out).exists(): raise ValueError('output exists')
 if m.get('code_sha256')!=sha(__file__) or m.get('python')!=sys.version.split()[0]: raise ValueError('code/runtime mismatch')
 if len(paths)!=6 or len(m.get('inputs',[]))!=6: raise ValueError('six source files required')
 for i,(v,r) in enumerate(zip(paths,m['inputs'])):
  f=Path(v); seconds=i*20; stamp=f'{seconds//60:02d}{seconds%60:02d}0'
  if path(v)!=r['path'] or r['satellite']!=19 or '_G19_s202632000'+stamp+'_e' not in f.name or sha(f)!=r['sha256']: raise ValueError('source identity mismatch')
 return m

def unpack_unsigned(dataset, bits):
 import numpy as np
 if dataset.dtype.kind!='i' or dataset.dtype.itemsize!=bits//8 or dataset.attrs.get('_Unsigned',b'') not in (b'true','true'): raise ValueError('unknown unsigned packing')
 raw=np.asarray(dataset[:],dtype=getattr(np,'int'+str(bits))); unsigned=raw.view(getattr(np,'uint'+str(bits))).astype(np.uint64)
 fill=dataset.attrs.get('_FillValue'); valid=np.ones(len(raw),dtype=bool)
 if fill is not None: valid &= raw!=int(np.asarray(fill).item())
 return unsigned,valid

def decode_file(f,index,file_sha):
 import numpy as np
 required=('group_time_offset','group_id','group_parent_flash_id','group_lat','group_lon','group_quality_flag')
 if any(k not in f for k in required): raise ValueError('missing required field')
 d=f['group_time_offset']; raw,valid_time=unpack_unsigned(d,16)
 unit=d.attrs.get('units',b'');unit=unit.decode() if isinstance(unit,bytes) else unit
 origin=dt.datetime.fromisoformat(unit.removeprefix('seconds since ')).replace(tzinfo=dt.timezone.utc) if unit.startswith('seconds since ') else None
 base=dt.datetime(2026,11,16,tzinfo=dt.timezone.utc)
 if origin is None or (origin-base).total_seconds()!=index*20: raise ValueError('time origin mismatch')
 scale=float(np.asarray(d.attrs['scale_factor']).item());off=float(np.asarray(d.attrs['add_offset']).item())
 if not math.isfinite(scale) or scale<=0 or not math.isfinite(off): raise ValueError('invalid packing')
 elapsed=raw*scale+off+index*20
 ids,vi=unpack_unsigned(f['group_id'],32); parents,vp=unpack_unsigned(f['group_parent_flash_id'],16); q,vq=unpack_unsigned(f['group_quality_flag'],16)
 lat=np.asarray(f['group_lat'][:],float);lon=np.asarray(f['group_lon'][:],float); n=len(raw)
 if any(len(x)!=n for x in (ids,parents,q,lat,lon)): raise ValueError('shape mismatch')
 rows=[]; counts={'attempted':n,'inside':0,'outside':0,'invalid_time':0,'good_qc':0,'degraded_qc':0,'unknown_or_fill_qc':0,'invalid_position':0,'invalid_id_or_parent':0}
 for i in range(n):
  t=float(elapsed[i]) if valid_time[i] and math.isfinite(elapsed[i]) else None
  inside=t is not None and 0<=t<120
  counts['invalid_time' if t is None else 'inside' if inside else 'outside']+=1
  qc=int(q[i]) if vq[i] else None
  counts['good_qc' if qc==0 else 'degraded_qc' if qc in (1,3,5) else 'unknown_or_fill_qc']+=1
  position=math.isfinite(lat[i]) and math.isfinite(lon[i]) and -90<=lat[i]<=90 and -180<=lon[i]<=180
  if not position: counts['invalid_position']+=1
  if not (vi[i] and vp[i]): counts['invalid_id_or_parent']+=1
  rows.append({'file_sha256':file_sha,'row_index':i,'elapsed_utc_s':t,'inside_interval':bool(inside),'group_id':int(ids[i]) if vi[i] else None,'parent_flash_id':int(parents[i]) if vp[i] else None,'lat_deg':float(lat[i]) if position else None,'lon_deg':float(lon[i]) if position else None,'quality':qc})
 return rows,counts

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-sha256',required=True);p.add_argument('--input',nargs=6,required=True);p.add_argument('--output',required=True);a=p.parse_args()
 m=authorize(a.manifest,a.authorization,a.manifest_sha256,a.authorization_sha256,a.input,a.output)
 try:
  import numpy,h5py
  if numpy.__version__!=m['numpy'] or h5py.__version__!=m['h5py']: raise ValueError('dependency version mismatch')
  out=[]
  for i,(p,row) in enumerate(zip(a.input,m['inputs'])):
   with h5py.File(p,'r') as f: rows,counts=decode_file(f,i,row['sha256'])
   out.append({'file':Path(p).name,'sha256':row['sha256'],'counts':counts,'rows':rows})
  result={'status':'completed','files':out,'limits':'Rows and IDs are file-local; duplicates retained. No duty cycle or physical event count established.'}
 except Exception as e: result={'status':'inconclusive','error_type':type(e).__name__,'reason':str(e)}
 result.update({'manifest_sha256':a.manifest_sha256,'authorization_sha256':a.authorization_sha256})
 with open(a.output,'x',encoding='utf-8',newline='\n') as f: json.dump(result,f,allow_nan=False,separators=(',',':'));f.write('\n')
 print(json.dumps({'status':result['status'],'file_counts':[x['counts'] for x in result.get('files',[])]}))
 if result['status']!='completed': raise SystemExit(2)
if __name__=='__main__':main()
