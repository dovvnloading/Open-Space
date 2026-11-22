"""Compare resumed synthetic multi-CCD searches to the untouched end-to-end search."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from astropy.io import fits
import numpy as np

import campaign
import scan_units as units
import search


def fixture(path):
    rng = np.random.default_rng(984)
    bg = 2000
    times = np.r_[rng.integers(0, 3000, bg), np.repeat(np.arange(1200,1300,4),2)].astype(float)
    chips = np.r_[rng.integers(0,2,bg), np.zeros(len(times)-bg,dtype=int)]
    cx, cy = rng.integers(3,1020,(2,len(times)))
    dx, dy = 10*np.sin(times/55), 12*np.cos(times/40)
    x, y = cx+dx, cy+dy
    x[bg:], y[bg:] = (512+rng.normal(0,.9,len(times)-bg) for _ in range(2))
    cx[bg:], cy[bg:] = np.rint(x[bg:]-dx[bg:]), np.rint(y[bg:]-dy[bg:])
    n = len(times)
    arr = [("TIME","D",times),("CHIPX","I",cx),("CHIPY","I",cy),("X","D",x),("Y","D",y),
           ("ENERGY","E",rng.uniform(900,2000,n)),("STATUS","32X",np.zeros((n,32),dtype=bool)),
           ("GRADE","I",np.zeros(n,dtype=int)),("CCD_ID","I",chips)]
    evt = fits.BinTableHDU.from_columns([fits.Column(name=k,format=f,array=v) for k,f,v in arr], name="EVENTS")
    evt.header.update({"INSTRUME":"ACIS","READMODE":"TIMED","GRATING":"NONE","OBS_ID":"98765",
                       "DTCOR":.96,"MJDREF":50814.,"TIMESYS":"TT","TIMEUNIT":"s","TIMEDEL":1.,
                       "TCTYP4":"RA---TAN","TCTYP5":"DEC--TAN","TCRPX4":512.,"TCRPX5":512.,
                       "TCRVL4":35.,"TCRVL5":-20.,"TCDLT4":-.492/3600,"TCDLT5":.492/3600,
                       "TCUNI4":"deg","TCUNI5":"deg"})
    gtis = []
    for chip, starts, stops in ((0,[0,2300],[2100,3000]),(1,[0],[3000])):
        gti = fits.BinTableHDU.from_columns([fits.Column(name="START",format="D",array=starts),
                                             fits.Column(name="STOP",format="D",array=stops)],name=f"GTI{chip}")
        gti.header["CCD_ID"] = chip
        gtis.append(gti)
    fits.HDUList([fits.PrimaryHDU(),evt,*gtis]).writeto(path)


class TestScanUnits(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        cls.root = Path(tmp.name)
        cls.path = cls.root / "acisf98765N001_evt2.fits"
        fixture(cls.path)
        evt = {"obsid":98765,"filename":cls.path.name,"sha256":units.file_sha(cls.path),
               "bytes":cls.path.stat().st_size,"url":"https://fixture.invalid/never-fetched"}
        obs = {"obsid":98765,"products":[{"input":evt,"complete":False}],"status":"error","acquired":[evt],
               "errors":["synthetic timeout"],"catalog_cache":"catalogs/fixture"}
        cls.ledger, cls.manifest = cls.root / "98765.json.gz", cls.root / "plan.json"
        campaign.save_obs(cls.ledger, obs)
        cls.receipt = units.plan(cls.ledger, cls.path, cls.manifest, 3)
        cls.baseline = search.search_file(cls.path, seed=cls.receipt["seed"], n_controls=24, n_injections=12)
        for i in reversed(range(3)):
            units.execute(cls.manifest, cls.path, cls.root / "parts" / f"group-{i}", group=i)

    def test_distributed_search_exactly_matches_stock_with_injection_indices(self):
        result = units.execute(self.manifest, self.path, self.root / "final", staged=self.root / "parts", ledger=self.ledger)
        self.assertEqual(units.encode(result), units.encode(self.baseline))
        self.assertTrue(result["triggers"])
        self.assertTrue(any(j["attempted"] for j in result["injections"]))
        self.assertEqual(self.receipt["attempted_injection_trials_in_call_order"],
                         [j["trial"] for j in self.baseline["injections"] if j["attempted"]])
        self.assertEqual(len(self.receipt["units"]), 3*25 + sum(j["attempted"] for j in self.baseline["injections"]))
        new = campaign.read_obs(self.root / "final/observations/98765.json.gz")
        self.assertEqual(new["products"][0]["search"], result)
        self.assertFalse(new["products"][0]["complete"])
        self.assertEqual(new["status"], "error")
        self.assertNotIn("search", campaign.read_obs(self.ledger)["products"][0])
        with self.assertRaisesRegex(ValueError, "already has a search"):
            units.product(new, self.path.name)

    def test_cached_group_skips_all_completed_scans(self):
        with patch.object(units, "ORIGINAL_SCAN", side_effect=AssertionError("completed unit ran again")):
            units.execute(self.manifest, self.path, self.root / "parts/group-1", group=1)

    def test_incomplete_or_disagreeing_artifacts_cannot_assemble(self):
        with self.assertRaisesRegex(ValueError, "Missing .* units"):
            units.execute(self.manifest, self.path, self.root / "partial", staged=self.root / "parts/group-0")
        one = next((self.root / "parts/group-0").rglob("*" + units.SUFFIX))
        with tempfile.TemporaryDirectory(dir=self.root) as d:
            folder = Path(d)
            for branch in ("a","b"):
                target = folder / branch / one.parent.name / one.name
                target.parent.mkdir(parents=True)
                shutil.copyfile(one,target)
            target.write_bytes(b"invalid duplicate")
            with self.assertRaisesRegex(ValueError, "disagree"):
                units.existing_files(folder,one.parent.name,{one.name.removesuffix(units.SUFFIX)})

    def test_changed_input_bytes_or_recorded_descriptor_cannot_be_used(self):
        with tempfile.TemporaryDirectory(dir=self.root) as d:
            path = Path(d)/self.path.name
            path.write_bytes(b"other file")
            with self.assertRaisesRegex(ValueError, "differs"):
                units.validate_plan(self.manifest,path)
        desc = self.receipt["units"][0]
        file, = (self.root / "parts").rglob(desc["id"]+units.SUFFIX)
        with self.assertRaisesRegex(ValueError, "corrupt"):
            units.payload(file, units.file_sha(self.manifest), {**desc,"rows":desc["rows"]+1})

    def test_atomic_retries_do_not_overwrite_different_results(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "record"
            with ThreadPoolExecutor(max_workers=4) as pool:
                list(pool.map(lambda _:units.atomic(path,b"first"), range(8)))
            with self.assertRaisesRegex(ValueError, "Conflicting"):
                units.atomic(path,b"second")
            self.assertEqual(path.read_bytes(),b"first")
