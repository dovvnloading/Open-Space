"""Test the original-ObsID matching rule and real/control audit symmetry."""
import unittest

from csc_variability import characterize, sky_predicate, triggers


class CSCAuditTests(unittest.TestCase):
    def test_real_and_control_have_equal_catalog_results(self):
        candidate={"ra_deg":0.001,"dec_deg":2.,"position_error_arcsec":1.2,"instrument_pass":True}
        other={**candidate,"instrument_pass":False}
        record={"products":[{"complete":True,"input":{"filename":"a","bytes":1,"sha256":"abc","url":"test"},
                             "search":{"triggers":[other,candidate],"controls":[{"replica":0,"triggers":[candidate]}]}}]}
        rows,_=triggers(record)
        self.assertEqual([r["trigger_key"] for r in rows],["real/1","control-0/0"])
        master={"source":{"name":"source","ra":"0.0011","dec":"2","err_ellipse_r0":"1"}}
        detection={"source":[{"match_type":"u","var_index_b":"10","flux_significance_b":"8.9"}]}
        a,b=[characterize(r,master,detection,"ok") for r in rows]
        self.assertEqual(a,b)
        self.assertTrue(a["matches"][0]["unique_association_var_index_at_least_8"])
        self.assertIn(" OR ",sky_predicate(rows))

    def test_failed_catalog_differs_from_empty_catalog(self):
        row={"ra_deg":210.5,"dec_deg":-45.,"position_error_arcsec":1.5}
        self.assertEqual(characterize(row,{},{},"incomplete")["status"],"incomplete_catalog_queries")
        self.assertEqual(characterize(row,{},{},"ok")["status"],"no_master_within_search_radius")


if __name__=="__main__":
    unittest.main()
