"""Guardrail: never use HTML error-page fragments as neutron station IDs."""
import importlib.util
from pathlib import Path
import unittest

p=Path(__file__).resolve().parents[1]/"examples"/"audit_neutron_monitor_identity.py"
s=importlib.util.spec_from_file_location("audit_neutron_monitor_identity",p)
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

class StationIdentityAuditTests(unittest.TestCase):
    def test_html_tokens_are_not_stations(self):
        for invalid in ('"-//W3C//DTD','4.01//EN"','<!DOCTYPE','HTML','PUBLIC'):
            self.assertTrue(m.looks_like_html_header(invalid),invalid)
        for valid in ('OULU','JUNG','SOPO','KIEL'):
            self.assertFalse(m.looks_like_html_header(valid),valid)

    def test_counts_and_caution_are_traceable(self):
        rows=[{'station':'"-//W3C//DTD','station_name':'"-//W3C//DTD','count_rate':100},
              {'station':'OULU','station_name':'Oulu, Finland','count_rate':None}]
        r=m.summarize(rows)
        self.assertEqual(r['row_count'],2)
        self.assertEqual(r['html_doctype_fragment_rows'],1)
        self.assertEqual(r['missing_count_rate_rows'],1)
        self.assertEqual(r['html_doctype_fragment_fraction'],.5)

if __name__=='__main__':
    unittest.main()
