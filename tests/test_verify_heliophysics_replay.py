"""Guard against silent upstream archive or statistic changes."""
import importlib.util
from pathlib import Path
import unittest

file = Path(__file__).resolve().parents[1] / "examples" / "verify_heliophysics_replay.py"
spec = importlib.util.spec_from_file_location("verify_heliophysics_replay", file)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class ArchiveReplayValidatorTests(unittest.TestCase):
    def test_missing_or_changed_fields_reported(self):
        self.assertEqual(mod.mismatches({"x":[1,2]}, {"x":[1,2]}), [])
        self.assertEqual(mod.mismatches({"x":1,"zero":0}, {"x":1.0,"zero":0.0}), [])
        self.assertTrue(mod.mismatches({"flag":True}, {"flag":1}))
        self.assertTrue(any("x[1]" in d for d in mod.mismatches({"x":[1,2]}, {"x":[1,3]})))
        self.assertTrue(any("missing" in d for d in mod.mismatches({"x":1}, {})))
        self.assertTrue(any("unexpected" in d for d in mod.mismatches({}, {"x":1})))

    def test_cme_source_lock_and_exact_output(self):
        report={"source_sha256":mod.SOURCE_SHA256["cme"],"cohorts":{"all":{"n":113,"r":-0.379398}}}
        mod.verify("cme", report, report)
        with self.assertRaisesRegex(ValueError,"Measurement packet differs"):
            mod.verify("cme", report, {**report,"cohorts":{"all":{"n":112,"r":-0.379398}}})
        with self.assertRaisesRegex(ValueError,"source digest"):
            mod.verify("cme", report, {**report,"source_sha256":"bad"})

    def test_forbush_lock_validates_both_independent_inputs(self):
        report={"source_files":{k:{"sha256":v} for k,v in mod.SOURCE_SHA256["forbush"].items()}}
        mod.verify("forbush", report, report)
        corrupted={"source_files":{"feid":{"sha256":mod.SOURCE_SHA256["forbush"]["feid"]},
                                   "omni":{"sha256":"wrong"}}}
        with self.assertRaisesRegex(ValueError,"source digest"):
            mod.verify("forbush", report, corrupted)

    def test_neutron_source_identity_lock(self):
        report={"input_sha256":mod.SOURCE_SHA256["neutron"],"html_doctype_fragment_rows":497835}
        mod.verify("neutron",report,report)
        with self.assertRaises(ValueError):
            mod.verify("neutron",report,{**report,"html_doctype_fragment_rows":0})


if __name__ == "__main__":
    unittest.main()
