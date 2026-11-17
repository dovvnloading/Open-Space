import hashlib,json,tempfile,unittest,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np
from astropy.io import fits
import development as d

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fixture(p,value,unit='nT'):
 h=fits.Header()
 for k,v in {'CTYPE1':'HPLN-TAN','CTYPE2':'HPLT-TAN','CUNIT1':'deg','CUNIT2':'deg','CRPIX1':1.,'CRPIX2':1.,'CRVAL1':0.,'CRVAL2':0.,'CDELT1':1.,'CDELT2':1.,'PC1_1':1.,'PC1_2':0.,'PC2_1':0.,'PC2_2':1.,'BUNIT':unit}.items(): h[k]=v
 i=fits.ImageHDU(np.full((3,3),value,dtype=float),header=h,name='PRIMARY DATA ARRAY')
 u=fits.ImageHDU(np.ones((3,3)),header=h,name='UNCERTAINTY ARRAY')
 fits.HDUList([fits.PrimaryHDU(),i,u]).writeto(p)
class Tests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns'); self.r=Path(self.tmp.name); self.x=self.r/'PUNCH_L2_CTM_20261101000029_v0l.fits'; self.y=self.r/'PUNCH_L2_CTM_20261101000829_v0l.fits'; fixture(self.x,1); fixture(self.y,3)
  self.m=self.r/'manifest.json'; self.m.write_text(json.dumps({'python':sys.version.split()[0],'code_sha256':sha(d.__file__),'output_path':str((self.r/'out.json').resolve()),'inputs':[{'name':p.name,'sha256':sha(p),'path':str(p.resolve())} for p in [self.x,self.y]]}))
  self.a=self.r/'auth.json'; self.a.write_text(json.dumps({'phase':'two-image-development','manifest_sha256':sha(self.m),'output_path':str((self.r/'out.json').resolve())})); self.out=self.r/'out.json'
 def tearDown(self): self.tmp.cleanup()
 def call(self): return d.authorize(self.m,self.a,sha(self.m),sha(self.a),[self.x,self.y],self.out)
 def test_actual_fits_difference(self):
  self.call(); result=d.characterize([self.x,self.y]); self.assertEqual(result['difference_median'],2); self.assertEqual(result['eligible_pixels'],9)
 def test_input_mutation_fails_before_decode(self):
  self.x.write_bytes(self.x.read_bytes()+b'x')
  with patch.object(d,'characterize',side_effect=AssertionError('decode')):
   with self.assertRaises(ValueError): self.call()
 def test_wrong_phase_and_output_collision(self):
  self.a.write_text(json.dumps({'phase':'holdout','manifest_sha256':sha(self.m),'output_path':str((self.r/'out.json').resolve())}))
  with self.assertRaises(ValueError): self.call()
  self.a.write_text(json.dumps({'phase':'two-image-development','manifest_sha256':sha(self.m),'output_path':str((self.r/'out.json').resolve())})); self.out.touch()
  with self.assertRaises(ValueError): self.call()
 def test_scope_and_auth_hash(self):
  with self.assertRaises(ValueError): d.authorize(self.m,self.a,sha(self.m),sha(self.a),[self.x,self.x],self.out)
  with self.assertRaises(ValueError): d.authorize(self.m,self.a,sha(self.m),'0'*64,[self.x,self.y],self.out)
 def test_wcs_mismatch_is_inconclusive(self):
  fixture(self.r/'third.fits',4,'other')
  with self.assertRaises(ValueError): d.characterize([self.x,self.r/'third.fits'])
 def test_wrong_path_and_authorized_output(self):
  copied=self.r/'other'/self.x.name; copied.parent.mkdir(); copied.write_bytes(self.x.read_bytes())
  with self.assertRaises(ValueError): d.authorize(self.m,self.a,sha(self.m),sha(self.a),[copied,self.y],self.out)
  with self.assertRaises(ValueError): d.authorize(self.m,self.a,sha(self.m),sha(self.a),[self.x,self.y],self.r/'other.json')
 def test_zero_uncertainty_excluded(self):
  with fits.open(self.x,mode='update') as h: h['UNCERTAINTY ARRAY'].data[0,0]=0
  result=d.characterize([self.x,self.y]); self.assertEqual(result['common_finite_pixels'],9); self.assertEqual(result['positive_uncertainty_pixels'],8); self.assertEqual(result['eligible_pixels'],8)
 def test_headers_checked_before_any_data(self):
  class NoData:
   def __init__(self,name,header): self.name=name; self.header=header
   @property
   def data(self): raise AssertionError('decoded before header gate')
  with fits.open(self.x) as h: a=NoData('PRIMARY DATA ARRAY',h['PRIMARY DATA ARRAY'].header.copy()); u=NoData('UNCERTAINTY ARRAY',h['UNCERTAINTY ARRAY'].header.copy())
  wrong=u.header.copy(); wrong['BUNIT']='wrong'
  with self.assertRaises(ValueError): d.inspect_headers([{'PRIMARY DATA ARRAY':a,'UNCERTAINTY ARRAY':u},{'PRIMARY DATA ARRAY':a,'UNCERTAINTY ARRAY':NoData('UNCERTAINTY ARRAY',wrong)}])
if __name__=='__main__': unittest.main()
