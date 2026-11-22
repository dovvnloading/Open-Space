"""Synthetic cases for the level-2 frame audit, run by the Chandra workflow."""
import csv
import gzip
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from astropy.io import fits
import diagnostics
import frame_forensics as ff


class FrameForensicsTests(unittest.TestCase):
    def table(self,has_frame=True):
        x=np.array([100.,100.1,100.2,100.3,100.4,100.5,100.2,110.,111.,180.,100.1])
        n=len(x)
        fld={"TIME":("D",[11.]*10+[13.]),"CCD_ID":("I",[7]*n),"X":("D",x),"Y":("D",[100.]*n),
             "CHIPX":("I",[100,100,100,101,101,101,101,110,111,180,102]),"CHIPY":("I",[500]*n),
             "GRADE":("I",[0]*n),"STATUS":("32X",np.zeros((n,32),dtype=bool)),"ENERGY":("E",[900.]*n)}
        fld["STATUS"][1][6,2]=True
        if has_frame:
            fld["EXPNO"]=("J",[10]*8+[11,10,12])
            fld["FLTGRADE"]=("I",[13]*n)
        return fits.BinTableHDU.from_columns([fits.Column(name=k,format=v[0],array=v[1]) for k,v in fld.items()])

    def test_exposure_number_dirty_rows_and_matched_controls(self):
        trig={"ccd_id":7,"aperture_x":100.,"aperture_y":100.,"aperture_radius_px":2.,
              "peak_start_met_s":10.,"peak_window_s":2.,"start_met_s":10.,"stop_met_s":15.}
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)
            m=ff.analyze(self.table(),trig,set(range(6)),8.,15.,out)
            self.assertEqual((m["peak_photons"],m["independent_exposure_numbers_in_peak"],m["observed_ccd_frame_keys"]),(6,1,3))
            self.assertEqual(m["peak_frames"][0]["key"],10)
            frame = m["node_row_patterns_across_trial"][0]
            self.assertEqual((frame["frame_key"],frame["source_peak_counts"]),(10,6))
            pattern, = frame["source_node_rows"]
            self.assertEqual((pattern["node"],pattern["chipy"],pattern["source_counts"]),(0,500,6))
            self.assertEqual(pattern["outside_aperture_level2_counts"],2)
            self.assertEqual(pattern["entire_row_chipx_range"],[100.,180.])
            self.assertEqual(m["peak_frames"][0]["control_counts"][0],1)
            self.assertEqual(m["peak_photon_geometry"]["max_on_one_chip_column"],3)
            self.assertEqual(m["peak_photon_geometry"]["max_on_one_chip_row"],6)
            self.assertEqual(m["source_occupancy_outside_peak_frames"],{"0":1,"1":1})
            self.assertEqual(m["control_frame_aperture_cells"],24)
            with gzip.open(out/"local_level2_rows.csv.gz","rt") as f:
                local={int(r["evt2_row_zero_based"]):r for r in csv.DictReader(f)}
            self.assertEqual(local[6]["survey_grade_energy_status_pass"],"0")
            self.assertEqual(local[6]["status_as_stored"].count("1"),1)
            self.assertEqual(local[7]["control_aperture"],"0")
            self.assertNotIn(9,local)
            with gzip.open(out/"full_peak_ccd_frames.csv.gz","rt") as f:
                peak={int(r["evt2_row_zero_based"]) for r in csv.DictReader(f)}
            self.assertIn(9,peak)
            self.assertNotIn(8,peak)
            m=ff.analyze(self.table(False),trig,set(range(6)),8.,15.,out)
            self.assertIsNone(m["independent_exposure_numbers_in_peak"])
            self.assertEqual(m["frame_key"],"identical_TIME_proxy")
            self.assertIn("FLTGRADE",m["absent_optional_columns"])
            with self.assertRaisesRegex(ValueError,"peak row IDs"):
                ff.analyze(self.table(),trig,{0},8.,15.,out)

    def test_original_event_rows_must_regenerate_exact_bytes(self):
        table=np.zeros(3,dtype=diagnostics.DTYPE)
        table["row"]=[9,50,81]
        table["t"]=[1,2,3]
        table["x"]=[0,4.5,0]
        table["e"]=600
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d)
            original=folder/"original.csv.gz"
            diagnostics._events_csv(original,table,np.array([True,False,True]),np.array([False,True,False]),1,4,1,2)
            prior={"aperture_x":0.,"aperture_y":0.,"aperture_radius_px":2.,"start_met_s":1,"stop_met_s":4,
                   "peak_start_met_s":1,"peak_window_s":1,"peak_counts":1,
                   "diagnostics":{"files":{"event_rows.csv.gz":ff.sha(original)}}}
            with patch.object(diagnostics,"_get_data",return_value=(table,0.,5.)):
                known,lo,hi=ff.verify_rows(Path("fake"),prior,original,folder)
                self.assertEqual((known,lo,hi),({9},0.,5.))
                original.write_bytes(b"different bytes")
                with self.assertRaisesRegex(ValueError,"byte for byte"):
                    ff.verify_rows(Path("fake"),prior,original,folder)
