"""Offline synthetic invariants for the Euclid DR1-F HI release parser and claims."""
import importlib.util
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("euclid_hi_run",HERE/"run.py")
pipeline=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=pipeline
spec.loader.exec_module(pipeline)

class FrozenEligibilityTests(unittest.TestCase):
    def test_alfalfa_real_column_names(self):
        x=pd.DataFrame([{"AGCNr":1234,"Name":"HI 1234","RAdeg_HI":12.5,
            "DECdeg_HI":-1.1,"Vhelio":3500,"W50":65,"HIflux":2.3,
            "Dist":47.1,"logMH":8.84,"HIcode":1}])
        out,rec=pipeline.normalize_catalog(x,{"id":"ALFALFA_alpha100","beam_fwhm_arcsec":210})
        self.assertEqual(out.iloc[0]["id"],"1234")
        self.assertAlmostEqual(out.iloc[0]["logMHI"],8.84)
        self.assertAlmostEqual(out.iloc[0]["flux_jykms"],2.3)
        self.assertEqual(rec["finite_logMHI"],1)
    def test_fashi_radio_columns(self):
        x=pd.DataFrame([{"name":"J123","fashi_id":11,"ra":13.1,"dec":4.2,"v_rad":3100,
             "W_50":50,"mass":9.2,"distance":42.3,"S_sum":2.1,"freq":1405.2}])
        out,rec=pipeline.normalize_catalog(x,{"id":"FASHI_DR2","beam_fwhm_arcsec":174})
        self.assertEqual(out.iloc[0]["id"],"11")
        self.assertAlmostEqual(out.iloc[0]["logMHI"],9.2)
        self.assertEqual(rec["finite_velocity"],1)
    def test_masked_missing_mass_not_zero(self):
        x=pd.DataFrame([{"AGCNr":1,"RAdeg_HI":25.,"DECdeg_HI":5.,
                         "Vhelio":2500,"W50":60,"HIflux":np.nan,"Dist":np.nan,"logMH":np.nan}])
        out,_=pipeline.normalize_catalog(x,{"id":"ALFALFA_alpha100","beam_fwhm_arcsec":210})
        self.assertTrue(np.isnan(out.iloc[0]["logMHI"]))
    def test_MOC_beam_edge_is_not_center(self):
        class LeftEdge:
            def contains_lonlat(self,lons,lats):
                return np.asarray(lons.to_value("deg"))>100
        a=pd.DataFrame([{"ra_deg":100.001,"dec_deg":0.,"beam_fwhm_arcsec":240.},
                        {"ra_deg":99.999,"dec_deg":0.,"beam_fwhm_arcsec":240.},
                        {"ra_deg":90.,"dec_deg":0.,"beam_fwhm_arcsec":240.}])
        r=pipeline.coverage(a,[LeftEdge()])
        self.assertEqual(r["euclid_coverage_code"].tolist(),[2,1,0])

if __name__=="__main__":
    unittest.main()
