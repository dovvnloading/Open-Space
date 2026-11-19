"""Small in-memory regression tests, independent of external publisher access."""
import math,unittest
import numpy as np
import run
class Tests(unittest.TestCase):
    def test_nominal_denominator(self):
        self.assertEqual(run.NP,1200)
        self.assertAlmostEqual(3*(run.NP*.1/3600)**2,.00333333333333,places=10)
        self.assertEqual({f["name"] for f in run.FROZEN["fields"]},{"north","south","fornax"})
    def test_robust(self):
        self.assertGreater(run.robust([1,2]*100),0)
        self.assertTrue(math.isnan(run.robust([1]*3)))
    def test_wcs(self):
        w=run.target_wcs(61.241,-48.423)
        x,y=w.all_world2pix([61.241],[-48.423],0)
        self.assertAlmostEqual(x[0],(run.NP-1)/2,places=3)
        self.assertAlmostEqual(y[0],(run.NP-1)/2,places=3)
    def test_url_validation(self):
        expected="dr1/MER/1/VIS/p.fits"
        self.assertEqual(run.key_from_url("https://irsacloud.ipac.caltech.edu/ibe/data/euclid/"+expected),expected)
        with self.assertRaises(ValueError):
            run.key_from_url("https://evil.invalid/ibe/data/euclid/"+expected)
    def test_local_psf_proxy_on_synthetic_sources(self):
        yy,xx=np.mgrid[0:run.NP,0:run.NP]
        a=np.zeros((run.NP,run.NP),dtype=np.float32)
        for x,y in [(200,200),(400,200),(600,200),(800,200),
                    (200,500),(400,500),(600,500),(800,500)]:
            a+=200*np.exp(-((xx-x)**2+(yy-y)**2)/(2*1.4**2))
        mask=np.ones_like(a,dtype=bool)
        rms=np.ones_like(a,dtype=np.float32)
        fit,rec=run.psf_calibrate(a,a*1.3,mask,rms,rms)
        self.assertIsNotNone(fit)
        self.assertEqual(rec["status"],"EMPIRICAL_COMPACT_SOURCE_PSF_AND_FLUX_PROXY")
        self.assertAlmostEqual(fit[2],1.3,places=3)
    def test_streak_rejected(self):
        a=np.zeros((run.NP,run.NP),dtype=np.float32)
        mask=np.zeros_like(a,dtype=bool)
        mask[200:1000,200:1000]=True
        a[600,250:950]=20
        found=run.peaks(a,mask,1,3)
        self.assertFalse(any(abs(q["x"]-600)<5 and abs(q["y"]-600)<5 for q in found))
if __name__=="__main__":
    unittest.main()
