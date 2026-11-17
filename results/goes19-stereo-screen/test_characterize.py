import json,hashlib,pathlib,tempfile,unittest
import numpy as np,h5py
import characterize as c
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns');self.root=pathlib.Path(self.tmp.name);self.f=self.root/'f.h5'
  with h5py.File(self.f,'w') as h:
   for k,v,dtype in [('group_time_offset',[0,10,20,30],np.int16),('group_id',[1,1,2,3],np.int32),('group_parent_flash_id',[4,4,5,6],np.int16),('group_quality_flag',[0,3,-1,7],np.int16)]:
    d=h.create_dataset(k,data=np.array(v,dtype=dtype));d.attrs['_Unsigned']='true'
    if k=='group_quality_flag': d.attrs['_FillValue']=np.int16(-1)
    if k=='group_time_offset': d.attrs['units']='seconds since 2026-11-16 00:00:00.000';d.attrs['scale_factor']=1.;d.attrs['add_offset']=-5.
   h.create_dataset('group_lat',data=[10.,10.,np.nan,10.]);h.create_dataset('group_lon',data=[-100.]*4)
 def tearDown(self):self.tmp.cleanup()
 def test_rows_time_flags_and_duplicates(self):
  with h5py.File(self.f) as h:rows,n=c.decode_file(h,0,'digest')
  self.assertEqual(n['attempted'],4);self.assertEqual(n['outside'],1);self.assertEqual(n['inside'],3);self.assertEqual(n['good_qc'],1);self.assertEqual(n['degraded_qc'],1);self.assertEqual(n['unknown_or_fill_qc'],2);self.assertEqual(n['invalid_position'],1);self.assertEqual(rows[0]['elapsed_utc_s'],-5);self.assertEqual([r['group_id'] for r in rows[:2]],[1,1]);self.assertEqual([r['row_index'] for r in rows],list(range(4)))
 def test_bad_origin(self):
  with h5py.File(self.f,'r+') as h:h['group_time_offset'].attrs['units']='seconds since 2026-11-16 00:00:20.000'
  with h5py.File(self.f) as h:
   with self.assertRaises(ValueError):c.decode_file(h,0,'x')
 def test_packing_and_time_fill(self):
  with h5py.File(self.f,'r+') as h:h['group_time_offset'].attrs['_FillValue']=np.int16(10)
  with h5py.File(self.f) as h:rows,n=c.decode_file(h,0,'x')
  self.assertEqual(n['invalid_time'],1);self.assertIsNone(rows[1]['elapsed_utc_s'])
 def test_phase_and_output_collision(self):
  import sys
  inputs=[];paths=[]
  for i in range(6):
   sec=i*20;stamp=f'{sec//60:02d}{sec%60:02d}0';p=self.root/('OR_GLM-L2-LCFA_G19_s202632000'+stamp+'_e.nc');p.write_bytes(self.f.read_bytes());paths.append(p);inputs.append({'satellite':19,'path':str(p.resolve()),'sha256':c.sha(p)})
  out=self.root/'out.json';m=self.root/'m.json';a=self.root/'a.json';m.write_text(json.dumps({'phase':'g19-characterization','python':sys.version.split()[0],'code_sha256':c.sha(c.__file__),'inputs':inputs,'output_path':str(out.resolve())}));a.write_text(json.dumps({'phase':'g19-characterization','manifest_sha256':c.sha(m),'output_path':str(out.resolve())}))
  c.authorize(m,a,c.sha(m),c.sha(a),paths,out);out.touch()
  with self.assertRaises(ValueError):c.authorize(m,a,c.sha(m),c.sha(a),paths,out)
  out.unlink()
  with self.assertRaises(ValueError):c.authorize(m,a,c.sha(m),c.sha(a),paths,self.root/'wrong.json')
  paths[0].write_bytes(b'changed')
  with self.assertRaises(ValueError):c.authorize(m,a,c.sha(m),c.sha(a),paths,out)
if __name__=='__main__':unittest.main()
