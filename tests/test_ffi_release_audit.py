"""Offline checks of the provisional LMC-East release and external-file auditor."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from zipfile import ZipFile

from open_space.ffi_release_audit import validate_report, audit_artifacts


class ImageReleaseAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[1]
        cls.report = json.loads((root / "results/ffi-lmc-east-preliminary.json").read_text())

    def test_documented_geometry_and_denominators(self):
        stats = validate_report(self.report)
        self.assertEqual(stats["planned_fields"], 6)
        self.assertEqual(stats["grid_centers_icrs_deg"], [
            [78.2, -69.4], [78.6, -69.4], [79.0, -69.4],
            [78.2, -69.2], [78.6, -69.2], [79.0, -69.2]])
        self.assertEqual(stats["recurrent_sources"], 13)
        self.assertEqual(stats["injection_recovery"]["1.0"]["fraction"], 13/18)
        self.assertAlmostEqual(stats["null_controls"]["half_day_block_shuffle"]["recurrences_per_sq_degree_trial"], 13/0.799)
        self.assertEqual(stats["scientific_reproduction"], "not_performed")

    def test_invalid_outcome_partition_is_rejected(self):
        altered = deepcopy(self.report)
        altered["reported_scan"]["first_sector_dispositions"]["single_sector"] += 1
        with self.assertRaisesRegex(ValueError, "do not sum"):
            validate_report(altered)

    def test_invalid_source_harmonic_is_rejected(self):
        altered = deepcopy(self.report)
        altered["reported_source_followup"][0]["harmonic_order"] = 2
        with self.assertRaisesRegex(ValueError, "harmonic"):
            validate_report(altered)

    def test_missing_artifacts_cannot_be_described_as_verified(self):
        with TemporaryDirectory() as tmp:
            check = audit_artifacts(self.report, tmp)
            self.assertEqual(check["state"], "external_evidence_incomplete")
            self.assertEqual(len(check["artifacts"]), 7)
            self.assertTrue(all(x["status"] == "not_available" for x in check["artifacts"]))

    def test_hashes_real_zip_and_checks_integrity(self):
        with TemporaryDirectory() as tmp:
            archive = Path(tmp) / "test.zip"
            with ZipFile(archive, "w") as z:
                z.writestr("observation.txt", "independent evidence\n")
            data = archive.read_bytes()
            altered = deepcopy(self.report)
            altered["external_artifacts"] = [{"filename": "test.zip", "bytes": len(data),
                                                "sha256": sha256(data).hexdigest()}]
            self.assertEqual(validate_report(altered)["external_artifacts_expected"], 1)
            result = audit_artifacts(altered, tmp, check_zip_crc=True)
            self.assertEqual(result["state"], "all_files_locally_hashed")
            self.assertTrue(result["artifacts"][0]["zip_crc_checked"])
            self.assertEqual(result["artifacts"][0]["sha256"], sha256(data).hexdigest())
            altered["external_artifacts"][0]["sha256"] = "0"*64
            self.assertEqual(audit_artifacts(altered, tmp)["artifacts"][0]["status"], "hash_mismatch")

    def test_rejects_unsafe_filenames_and_bool_counts(self):
        altered = deepcopy(self.report)
        altered["external_artifacts"][0]["filename"] = "../danger.zip"
        with self.assertRaisesRegex(ValueError, "safe basenames"):
            validate_report(altered)
        altered = deepcopy(self.report)
        altered["reported_scan"]["first_sector_peaks"] = True
        with self.assertRaisesRegex(ValueError, "integer"):
            validate_report(altered)


if __name__ == "__main__":
    unittest.main()
