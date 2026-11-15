"""Offline frozen-table audit and independent-paper diagnostic regression tests."""
import csv
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / 'results' / 'oec-orbital-scale-audit'
spec = importlib.util.spec_from_file_location('oec_audit', HERE / 'audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def example_table():
    # Source uses quoted multi-list entries, which must not split into columns.
    return io.StringIO(
        'name,binaryflag,period,semimajoraxis,hoststar_mass,hoststar_radius,discoverymethod,lastupdate,list\n'
        '51 Peg b,0,4.231,0.052,1.04,1.2,RV,15/05/24,Confirmed planets\n'
        'Kepler-445 b,0,2.984151,0.102,0.18,0.21,transit,15/01/06,Confirmed planets\n'
        'Kepler-445 c,0,4.871229,0.141,0.18,0.21,transit,15/01/06,"Confirmed planets, Multiples"\n'
        'Kepler-441 b,0,207.2482,0.64,0.0572,0.55,transit,15/01/28,Confirmed planets\n'
        'not confirmed,0,3,0.03,0.9,1.0,transit,15/01/01,Controversial\n'
        'not single,2,3,0.03,0.9,1.0,transit,15/01/01,Confirmed planets\n'
        'missing axis,0,3,,0.9,1.0,transit,15/01/01,Confirmed planets\n'
        'invalid infinity,0,3,inf,0.9,1.0,transit,15/01/01,Confirmed planets\n'
    )


class OrbitAuditTests(unittest.TestCase):
    def test_fixed_analytic_controls(self):
        self.assertAlmostEqual(audit.implied_mass_solar(365.25, 1), 1)
        self.assertAlmostEqual(audit.predicted_axis_au(365.25, 1), 1)
        self.assertAlmostEqual(audit.SOLAR_RADIUS_IN_AU, 0.004650467, places=7)
        self.assertAlmostEqual(audit.implied_mass_solar(365.25, 2), 8)
        self.assertAlmostEqual(audit.predicted_axis_au(365.25, 8), 2)

    def test_audit_and_exclusions(self):
        flags, counts = audit.audit_stream(example_table())
        self.assertEqual(counts['source_rows'], 8)
        self.assertEqual(counts['eligible'], 4)
        self.assertEqual(counts['excluded_not_confirmed'], 1)
        self.assertEqual(counts['excluded_nonzero_binaryflag'], 1)
        self.assertEqual(counts['excluded_incomplete_or_invalid'], 2)
        self.assertEqual(len(flags), 3)
        self.assertEqual(flags[0]['name'], 'Kepler-445 b')
        self.assertEqual(set(f['name'] for f in flags),
                         {'Kepler-445 b', 'Kepler-445 c', 'Kepler-441 b'})
        self.assertAlmostEqual(float(flags[0]['mass_ratio_q']), 88.3216, places=3)

    def test_threshold_and_non_finite(self):
        with self.assertRaises(ValueError):
            audit.audit_stream(example_table(), 1)
        with self.assertRaises(ValueError):
            audit.audit_stream(example_table(), float('nan'))
        with self.assertRaises(ValueError):
            audit.audit_stream(io.StringIO('name,period\nA,1\n'))
        flags, counts = audit.audit_stream(example_table(), 100)
        self.assertEqual(flags, [])
        self.assertEqual(counts['flags_over_threshold'], 0)

    def test_primary_paper_mechanism_and_independent_irradiation(self):
        with (HERE / 'primary_checks.csv').open(newline='', encoding='utf-8') as f:
            checks = audit.primary_paper_diagnostics(f)
        self.assertEqual(len(checks), 6)
        for check in checks:
            self.assertLess(abs(float(check['oec_to_paper_scale_error_percent'])), 1.6)
            self.assertLess(abs(float(check['paper_flux_error_percent'])), 1.2)
            self.assertLess(float(check['using_oec_a_flux_earth']),
                            float(check['paper_flux_earth']) / 16)
            self.assertGreater(float(check['apparent_mass_ratio_q']), 70)

    def test_k2_285_source_paper_not_oec_transcription(self):
        with (HERE / 'k2_285_primary_checks.csv').open(newline='', encoding='utf-8') as f:
            cases = audit.source_paper_consistency(f)
        self.assertEqual(len(cases), 4)
        expected_m = [0.6155, 1.4649, 1.9948, 3.5943]
        for case, target_mass in zip(cases, expected_m):
            self.assertAlmostEqual(float(case['mass_implied_by_p_and_a_solar']),
                                   target_mass, delta=0.001)
            self.assertLess(abs(float(case['paper_axis_to_ratio_axis_percent'])), 0.5)
        self.assertLess(float(cases[0]['mass_ratio_to_paper_host']), 0.75)
        self.assertGreater(float(cases[3]['mass_ratio_to_paper_host']), 4.3)

    def test_cli_output_sha_and_deterministic_serialization(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            source = work / 'catalog.csv'
            source.write_text(example_table().getvalue(), encoding='utf-8')
            out = work / 'report'
            self.assertEqual(audit.main(['--input', str(source), '--out', str(out),
                                         '--primary-checks', str(HERE / 'primary_checks.csv'),
                                         '--source-model-checks', str(HERE / 'k2_285_primary_checks.csv')]), 0)
            summary = json.loads((out / 'summary.json').read_text())
            self.assertEqual(summary['counts']['source_rows'], 8)
            self.assertEqual(summary['independently_checked_primary_paper_rows'], 6)
            self.assertEqual(summary['source_paper_model_rows_checked'], 4)
            with (out / 'flags.csv').open(newline='', encoding='utf-8') as f:
                self.assertEqual(len(list(csv.DictReader(f))), 3)
            with (out / 'unit_conversion_checks.csv').open(newline='', encoding='utf-8') as f:
                self.assertEqual(len(list(csv.DictReader(f))), 6)
            # Reprocessing identical inputs yields identical CSV bytes.
            previous = (out / 'flags.csv').read_bytes()
            audit.main(['--input', str(source), '--out', str(out)])
            self.assertEqual((out / 'flags.csv').read_bytes(), previous)

if __name__ == '__main__':
    unittest.main()
