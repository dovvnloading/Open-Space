from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from open_space import tess, tess_archive


def listing(tic, sector):
    name = f"tess2018001000000-s{sector:04d}-{tic:016d}-0120-s_lc.fits"
    return f"curl -C - -L -o {name} https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/{name}\n"


class TessArchiveTests(unittest.TestCase):
    def test_inventory_intersection_pages(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sector, tics in ((1, [1, 9, 3, 7]), (2, [9, 3, 7, 8])):
                (root/f"s{sector}.sh").write_text("malformed input\n" + "".join(listing(t,sector) for t in tics))
            one = tess_archive.make_plan(1,2,2,0,"sample",root)
            two = tess_archive.make_plan(1,2,2,2,"sample",root)
            ids = [t["tic_id"] for t in one["targets"] + two["targets"]]
            expected = sorted([3,7,9], key=lambda t: sha256(f"sample:{t}".encode()).hexdigest())
            self.assertEqual(ids,expected)
            self.assertEqual(one["selection"]["shared_tics"],3)
            with self.assertRaises(ValueError):
                tess_archive.make_plan(1,2,2,4,"sample",root)

    def test_archived_sample_accounts_for_all_targets(self):
        root = Path(__file__).resolve().parents[1]
        plan = root/"data/tess-plan.json"
        result = json.loads((root/"results/tess-pilot.json").read_bytes())
        self.assertEqual(result["plan_sha256"],sha256(plan.read_bytes()).hexdigest())
        target = result["targets"][0]
        self.assertEqual(target["tic_id"],100100827)
        self.assertLess(abs(target["ephemeris"]["period_days"]/0.9414522-1),0.001)
        self.assertEqual(len(result["targets"]),25)
        self.assertEqual(result["dispositions"]["periodic_dip_with_veto"],2)

    @unittest.skipUnless(importlib.util.find_spec("astropy"), "archive extra required")
    def test_corrupt_cache_fails_before_analysis(self):
        with TemporaryDirectory() as tmp:
            name = listing(7,1).split()[5]
            p = {"filename":name,"sector":1,"uri":"mast:TESS/product/"+name,"sha256":"0"*64}
            (Path(tmp)/name).write_bytes(b"damaged")
            with self.assertRaisesRegex(ValueError,"Checksum"):
                tess_archive.cached_product(p,Path(tmp))


@unittest.skipUnless(importlib.util.find_spec("astropy"), "archive extra required")
class TessScienceTests(unittest.TestCase):
    @staticmethod
    def mock_curve():
        import numpy as np
        t = np.arange(0,18,2/1440)
        t = t[(t<8.2)|(t>9.7)]
        rng = np.random.default_rng(904)
        flux = 10000*(1+0.004*np.cos(2*np.pi*t/14))
        flux *= 1 - 0.013 * (abs((t-0.44+1.25)%2.5-1.25)<1/24)
        flux += rng.normal(size=len(t))*6
        return t,flux,np.ones(len(t))*6

    def test_known_period_and_missing_sector_preserves_gaps(self):
        import numpy as np
        curve=tess.flatten(self.mock_curve())
        self.assertTrue(np.all((curve[0]<8.2)|(curve[0]>9.7)))
        e=tess.blind(curve)
        self.assertAlmostEqual(e["period_days"],2.5,delta=0.01)
        fixed=tess.assess(curve,e)
        self.assertGreaterEqual(fixed["supported_events"],5)
        empty=tuple(np.array([]) for _ in range(3))
        self.assertEqual(tess.search_pair(curve,empty)["disposition"],"insufficient_coverage")

    def test_local_estimator_rejects_uncovered_prediction(self):
        import numpy as np
        t=np.arange(100,102,1/144)
        f=np.zeros(len(t))
        e=np.ones(len(t))*0.001
        ephem={"period_days":4,"duration_days":1/24,"epoch_btjd":103}
        result=tess.assess((t,f,e),ephem)
        self.assertEqual(result["covered_events"],0)
        self.assertIsNone(result["depth"])


@unittest.skipUnless(importlib.util.find_spec("astropy"), "archive extra required")
class MorphologyTests(unittest.TestCase):
    def test_heldout_box_beats_smooth_shapes_for_clean_transits(self):
        import numpy as np
        from open_space.tess_morphology import evaluate

        epoch, period, duration = 1.0, 2.0, 0.12
        def curve(start):
            t = np.arange(start, start + 18, 0.01)
            inside = abs((t - epoch + period/2) % period - period/2) <= duration/2
            return t, -0.02 * inside.astype(float), np.full(len(t), 0.001)
        ephem = {"epoch_btjd": epoch, "period_days": period, "duration_days": duration}
        result = evaluate(curve(0), curve(22), ephem)
        self.assertEqual(result["best_test_shape"], "box")
        self.assertEqual(result["test_samples"], 1800)


class CounterpartAndLaterTests(unittest.TestCase):
    def test_catalog_wraparound_and_exact_tic(self):
        from open_space.tess_eb_catalog import lookup
        rows = [{"tess_id": "1", "ra": "359.999", "dec": "0", "period": "2.9"},
                {"tess_id": "2", "ra": "5", "dec": "0", "period": "1.1"}]
        matched = lookup(rows, 2, 0.001, 0, 10)
        self.assertEqual([v["catalog_tic"] for v in matched], [1, 2])
        self.assertAlmostEqual(matched[0]["separation_arcsec"], 7.2, places=5)
        self.assertTrue(matched[1]["same_id"])

    @unittest.skipUnless(importlib.util.find_spec("astropy"), "archive extra required")
    def test_later_sinusoid_and_missing_product(self):
        import numpy as np
        from open_space.tess_reobserve import run, sinusoid
        t = np.arange(10, 38, 0.02)
        result = sinusoid((t, np.cos(2*np.pi*t/0.503), np.full(t.size, 0.02)))
        self.assertAlmostEqual(result["period_days"], 0.503, places=4)
        self.assertGreater(result["standard_power"], 0.99)
        record = run({"targets": [{"tic_id": 1, "product": None}]}, [], ".")
        self.assertEqual(record["cases"][0]["outcome"], "no_later_product")
