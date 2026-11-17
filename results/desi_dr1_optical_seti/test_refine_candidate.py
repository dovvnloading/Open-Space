import unittest
from laser_search import inject_line,scan_exposure
from refine_candidate import refit
from test_laser_search import fixture


class RefitTests(unittest.TestCase):
    def test_same_window_preserves_exact_subpixel_shapes(self):
        e=fixture()
        for phase in (.053,.51,.945,.955,.975,.995):
            true=e.wave[80]+phase*(e.wave[81]-e.wave[80])
            planted=inject_line(e,true,2000.)
            hits,_=scan_exposure(planted)
            v1=max(hits,key=lambda h:h.snr)
            result=refit(planted,v1.observed_A)
            self.assertAlmostEqual(result['bary_A'],true,places=7)
            self.assertAlmostEqual(result['flux17'],2000.,places=5)
            self.assertGreater(result['shape_p'],.999999)
