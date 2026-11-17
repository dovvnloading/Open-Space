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
  self.tmp=tempfile.TemporaryDirectory(); self.r=Path(self.tmp.name); self.x=self.r/'PUNCH_L2_CTM_20261101000029_v0l.fits'; self.y=self.r/'PUNCH_L2_CTM_20261101000829_v0l.fits'; fixture(self.x,1); fixture(self.y,3)
  self.m=self.r/'manifest.json'; self.m.write_text(json.dumps({'python':sys.version.split()[0],'code_sha256':sha(d.__file__),'inputs':[{'name':p.name,'sha256':sha(p)} for p in [self.x,self.y]]}))
  self.a=self.r/'auth.json'; self.a.write_text(json.dumps({'phase':'two-image-development','manifest_sha256':sha(self.m)})); self.out=self.r/'out.json'
 def tearDown(self): self.tmp.cleanup()
 def call(self): return d.authorize(self.m,self.a,sha(self.m),sha(self.a),[self.x,self.y],self.out)
 def test_actual_fits_difference(self):
  self.call(); result=d.characterize([self.x,self.y]); self.assertEqual(result['difference_median'],2); self.assertEqual(result['common_finite_pixels'],9)
 def test_input_mutation_fails_before_decode(self):
  self.x.write_bytes(self.x.read_bytes()+b'x')
  with patch.object(d,'characterize',side_effect=AssertionError('decode')):
   with self.assertRaises(ValueError): self.call()
 def test_wrong_phase_and_output_collision(self):
  self.a.write_text(json.dumps({'phase':'holdout','manifest_sha256':sha(self.m)}))
  with self.assertRaises(ValueError): self.call()
  self.a.write_text(json.dumps({'phase':'two-image-development','manifest_sha256':sha(self.m)})); self.out.touch()
  with self.assertRaises(ValueError): self.call()
 def test_scope_and_auth_hash(self):
  with self.assertRaises(ValueError): d.authorize(self.m,self.a,sha(self.m),sha(self.a),[self.x,self.x],self.out)
  with self.assertRaises(ValueError): d.authorize(self.m,self.a,sha(self.m),'0'*64,[self.x,self.y],self.out)
 def test_wcs_mismatch_is_inconclusive(self):
  fixture(self.r/'third.fits',4,'other')
  with self.assertRaises(ValueError): d.characterize([self.x,self.r/'third.fits'])
if __name__=='__main__': unittest.main()
