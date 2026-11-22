"""Synthetic EVT2 regression, run in the campaign Actions workflow."""
import json
from pathlib import Path
import tempfile
import unittest

from astropy.io import fits
import numpy as np

from diagnostics import save_diagnostics
from search import _time_support, search_file


class BlindSearchTests(unittest.TestCase):
    def test_union_with_overlapping_chip_livetime(self):
        self.assertEqual(_time_support([(0,5,.8),(2,7,.9),(8,9,.6)]),
                         [[0,2,.8],[2,7,.9],[8,9,.6]])

    def test_fits_to_trigger_with_controls_and_veto(self):
        rng = np.random.default_rng(49)
        # Two flares in an otherwise source-free, dithered field. The second
        # matches a single detector pixel; sky and detector tests must differ.
        nb, ns = 10000, 50
        time = np.r_[rng.integers(0,18000,nb), np.arange(5000,5150,3),
                     np.arange(9000,9150,3)].astype(float)
        dx = 8*np.sin(time/150)
        dy = 8*np.sin(time/125)
        cx, cy = rng.integers(2,1023,(2,nb+2*ns))
        x, y = cx + dx, cy + dy
        src, hot = slice(nb,nb+ns), slice(nb+ns,None)
        x[src], y[src] = 512+rng.normal(0,.75,ns), 512+rng.normal(0,.75,ns)
        cx[src], cy[src] = np.rint(x[src]-dx[src]), np.rint(y[src]-dy[src])
        cx[hot], cy[hot] = 750, 750
        x[hot], y[hot] = 750+dx[hot], 750+dy[hot]
        n = len(time)
        values = [("time","D",time),("chipx","I",cx),("chipy","I",cy),
                  ("x","D",x),("y","D",y),("energy","E",np.full(n,1500.)),
                  ("status","32X",np.zeros((n,32),dtype=bool)),
                  ("grade","I",np.zeros(n,dtype=int)),("ccd_id","I",np.zeros(n,dtype=int))]
        evt = fits.BinTableHDU.from_columns(
            [fits.Column(name=k,format=f,array=v) for k,f,v in values],name="EVENTS")
        evt.header.update({"INSTRUME":"ACIS","GRATING":"NONE","READMODE":"TIMED",
                           "OBS_ID":"99999","DTCOR":1.0,"MJDREF":50814.,"TIMESYS":"TT",
                           "TIMEUNIT":"s","TIMEDEL":1.,
                           "TCTYP4":"RA---TAN","TCTYP5":"DEC--TAN",
                           "TCRPX4":512.,"TCRPX5":512.,"TCRVL4":35.,"TCRVL5":-20.,
                           "TCDLT4":-.492/3600,"TCDLT5":.492/3600,
                           "TCUNI4":"deg","TCUNI5":"deg"})
        gti = fits.BinTableHDU.from_columns([
            fits.Column(name="START",format="D",array=[0.]),
            fits.Column(name="STOP",format="D",array=[18000.])],name="GTI")
        gti.header["CCD_ID"]=0
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"acisf99999_evt2.fits"
            fits.HDUList([fits.PrimaryHDU(),evt,gti]).writeto(path)
            result = search_file(path,seed=314,n_controls=3,n_injections=3)
            survivor = next(t for t in result["triggers"] if t["instrument_pass"])
            receipt = save_diagnostics(path, survivor, Path(directory)/"review")
            self.assertEqual(receipt["counts"]["on_counts"],survivor["on_counts"])
            self.assertEqual(set(receipt["files"]),{"diagnostics.png","lightcurves.csv","event_rows.csv.gz"})
        json.dumps(result,allow_nan=False)
        self.assertEqual(result["start_support_by_duration"]["60"], [[60,17880,1.]])
        self.assertEqual(result["exposure_s"],17820)
        real = [t for t in result["triggers"] if t["instrument_pass"]]
        self.assertTrue(any(abs(t["ra_deg"]-35)<.001 and t["before_s"]>=60 for t in real))
        self.assertTrue(any("repeated_detector_pixel" in t["vetoes"] for t in result["triggers"]))
        self.assertFalse(any(t["instrument_pass"] for c in result["controls"] for t in c["triggers"]))
        self.assertTrue(all(row["attempted"] for row in result["injections"]))


if __name__ == "__main__":
    unittest.main()
