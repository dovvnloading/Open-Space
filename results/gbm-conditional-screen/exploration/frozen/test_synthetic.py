import unittest,tempfile,sys,json,hashlib
from unittest.mock import patch
from pathlib import Path
import numpy as np
from astropy.io import fits
from simulator import simulate,score,monte_carlo,arrival_counts,bin_modulation,DT
from fits_reader import validate
import observational,manifest,exploratory
class Synthetic(unittest.TestCase):
 def test_joint_shared_latent_and_injection(self):
  a=monte_carlo(11,12);b=monte_carlo(12,12,injection=.25);self.assertGreater(np.median(b),np.median(a))
 def test_deadtime_reduces_counts(self):
  x=np.ones(256)*100000;a=arrival_counts(np.random.default_rng(3),x,0);b=arrival_counts(np.random.default_rng(3),x,2.6e-6);self.assertLess(b.sum(),a.sum());c=arrival_counts(np.random.default_rng(3),x,2.6e-6,overflow_prob=1);self.assertLess(c.sum(),b.sum())
 def test_envelope_alternative_finite(self):
  a,b=simulate(np.random.default_rng(1),kind='pulse_train');self.assertTrue(np.isfinite(score(a,b)))
 def test_poisson_normalization_nyquist(self):
  rng=np.random.default_rng(77);x=rng.poisson(20,size=(4000,256))-20;p=2*np.abs(np.fft.rfft(x,axis=1))**2/(256*20);self.assertAlmostEqual(float(p[:,40].mean()),2.,delta=.1);self.assertAlmostEqual(float((p[:,-1]*.5).mean()),1.,delta=.08)
 def test_bin_integrated_nyquist(self):
  t=np.arange(4)*DT;a=bin_modulation(t,32,0);np.testing.assert_allclose(a,[2/np.pi,-2/np.pi,2/np.pi,-2/np.pi],atol=1e-12);self.assertAlmostEqual(float(abs(bin_modulation(t,32,np.pi/2)).max()),0.,places=12)
 def test_fits_schema_energy_and_time(self):
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as tmp:
   p=Path(tmp)/'synthetic.fit';primary=fits.PrimaryHDU();primary.header['TIMESYS']='TT';primary.header['TIMEUNIT']='s';primary.header['DETNAM']='NAI_00';primary.header['OBJECT']='GRB000000000';primary.header['TRIGTIME']=100.;primary.header['MJDREFI']=51910;primary.header['MJDREFF']=.0007428703703703703;events=fits.BinTableHDU.from_columns([fits.Column(name='TIME',format='D',unit='s',array=[100.,101.,102.]),fits.Column(name='PHA',format='I',array=[0,1,2])],name='EVENTS');eb=fits.BinTableHDU.from_columns([fits.Column(name='CHANNEL',format='I',array=[0,1,2]),fits.Column(name='E_MIN',format='E',unit='keV',array=[10.,50.,100.]),fits.Column(name='E_MAX',format='E',unit='keV',array=[50.,100.,300.])],name='EBOUNDS');gti=fits.BinTableHDU.from_columns([fits.Column(name='START',format='D',unit='s',array=[99.]),fits.Column(name='STOP',format='D',unit='s',array=[104.])],name='GTI');fits.HDUList([primary,eb,events,gti]).writeto(p);t,_,_=validate(p,100,0,3);np.testing.assert_equal(t,[1.,2.]);
   with self.assertRaisesRegex(ValueError,'covered'):validate(p,100,0,5)
   with self.assertRaisesRegex(ValueError,'identity'):validate(p,100,0,3,detector='NAI_01')
   with self.assertRaisesRegex(ValueError,'identity'):validate(p,100,0,3,trigger_id='123456789')
   with fits.open(p,mode='update') as h:h['EBOUNDS'].columns['E_MIN'].unit='eV';h.flush()
   with self.assertRaisesRegex(ValueError,'energy units'):validate(p,100,0,3)
 def test_guard_mutation_and_collision(self):
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as tmp:
   p=Path(tmp)/'manifest.json';p.write_text('{}');digest=hashlib.sha256(p.read_bytes()).hexdigest();p.write_text('{"mutated":true}')
   with patch.object(manifest,'snapshot',return_value={}):
    with self.assertRaisesRegex(ValueError,'manifest mismatch'):manifest.verify(p,digest)
   out=Path(tmp)/'existing.json';out.write_text('existing');args=['x','--manifest','x','--manifest-sha256','x','--authorization','x','--authorization-sha256','x','--output',str(out)]
   with patch.object(sys,'argv',args),patch.object(observational,'verify',return_value={}),patch.object(observational,'authorize'),patch.object(observational.fits,'open',side_effect=AssertionError('opened')):
    with self.assertRaisesRegex(ValueError,'collision'):observational.main()
 def test_missing_auth_before_fits_open(self):
  args=['x','--manifest','x','--manifest-sha256','x','--authorization','x','--authorization-sha256','x','--output','D:/unused']
  with patch.object(sys,'argv',args),patch.object(observational,'verify',return_value={}),patch.object(observational,'authorize',side_effect=ValueError('sealed')),patch.object(observational.fits,'open',side_effect=AssertionError('decoded')):
   with self.assertRaisesRegex(ValueError,'sealed'):observational.main()
 def test_input_file_mutation(self):
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as tmp:
   root=Path(tmp);raw=root/'raw';raw.mkdir();rows=[]
   for i in range(4):
    name=f'file{i}.fit';(raw/name).write_bytes(b'original');rows.append({'url':'https://example.test/'+name,'sha256':hashlib.sha256(b'original').hexdigest()})
   receipt=root/'receipts.json';receipt.write_text(json.dumps({'files':rows}))
   def git(*args,**kwargs):return 'HEAD' if kwargs.get('text') and 'rev-parse' in args[0] else ''
   with patch.object(manifest,'RAW',raw),patch.object(manifest,'RECEIPTS',receipt),patch.object(manifest.subprocess,'check_output',side_effect=git):
    first=manifest.snapshot();self.assertEqual(len(first['tte_sha256']),4);(raw/'file1.fit').write_bytes(b'mutated')
    with self.assertRaisesRegex(ValueError,'input differs'):manifest.snapshot()
 def test_exploratory_fullband_deadtime_first(self):
  t=np.array([.1,.2,.3,1.1]);pha=np.array([1,127,2,127]);r=exploratory.describe(t,pha,np.array([1,2]),0.,1);self.assertEqual(r['fullband_count'],4);self.assertEqual(r['overflow_count'],2);self.assertEqual(r['windows'][0]['band_count'],2);self.assertAlmostEqual(r['windows'][0]['nominal_deadtime_seconds'],2*2.6e-6+2*1e-5);self.assertEqual(len(r['fourier_power'][0]),109)
 def test_exploratory_phase_guard(self):
  args=['x','--manifest','x','--manifest-sha256','x','--authorization','x','--authorization-sha256','x','--output','D:/unused']
  with patch.object(sys,'argv',args),patch.object(exploratory,'verify',return_value={}),patch.object(exploratory,'authorize',side_effect=ValueError('sealed')),patch.object(exploratory.fits,'open',side_effect=AssertionError('opened')):
   with self.assertRaisesRegex(ValueError,'sealed'):exploratory.main()
 def test_envelope_nuisance_and_refit(self):
  from envelope import fit,joint_bootstrap
  rng=np.random.default_rng(123);x=np.stack([rng.poisson(30,512),rng.poisson(20,512)]);y=np.stack([rng.poisson(25,1280),rng.poisson(18,1280)]);self.assertEqual(fit(x,33)['expected'].shape,x.shape)
  z=joint_bootstrap(np.random.default_rng(4),[x,y],3,33);self.assertEqual(len(z['maxima']),3);self.assertGreaterEqual(z['p'],.25)
if __name__=='__main__':unittest.main()
