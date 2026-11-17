import unittest,tempfile,pathlib,json,hashlib
import numpy as np,h5py
import reader
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns');self.root=pathlib.Path(self.tmp.name);self.paths=[]
  for i in range(6):
   sec=i*20; stamp=f'{sec//60:02d}{sec%60:02d}0'; p=self.root/('synthetic_G19_s202632000'+stamp+'_e.nc');self.paths.append(p)
   with h5py.File(p,'w') as h:
    x=h.create_dataset('group_time_offset',data=np.array([-1,0],dtype=np.int16));x.attrs['_Unsigned']='true';x.attrs['units']='seconds since 2026-11-16 00:00:00.000';x.attrs['scale_factor']=0.00038148;x.attrs['add_offset']=-5.
    for k,v in {'group_lat':[10.,11.],'group_lon':[-100.,-101.],'group_quality_flag':[0,3],'group_energy':[1,1],'group_parent_flash_id':[1,2]}.items(): h.create_dataset(k,data=v)
  sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
  self.m=self.root/'manifest.json';self.m.write_text(json.dumps({'phase':'source-development','reader_sha256':sha(pathlib.Path(reader.__file__)),'inputs':[{'satellite':19,'path':str(p.resolve()),'sha256':sha(p)} for p in self.paths]}));self.md=sha(self.m)
  self.a=self.root/'auth.json';self.a.write_text(json.dumps({'phase':'source-development','manifest_sha256':self.md}));self.ad=sha(self.a)
 def tearDown(self): self.tmp.cleanup()
 def test_real_hdf_fixture_unsigned(self):
  x=reader.read_source(self.m,self.a,self.md,self.ad,self.paths);self.assertAlmostEqual(x[0]['time_offset'][0],65535*0.00038148-5)
 def test_mutation_and_holdout_denied(self):
  self.paths[0].write_bytes(self.paths[0].read_bytes()+b'x')
  with self.assertRaises(ValueError):reader.read_source(self.m,self.a,self.md,self.ad,self.paths)
 def test_schema_units_rejected(self):
  with h5py.File(self.paths[0],'r+') as h: h['group_time_offset'].attrs['units']='milliseconds'
  with h5py.File(self.paths[0]) as h:
   with self.assertRaises(ValueError):reader.decode_groups(h)
if __name__=='__main__':unittest.main()
