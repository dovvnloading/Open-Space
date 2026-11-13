import copy
import json
from pathlib import Path
import tempfile
import unittest

from open_space.recurrence import packet_files, summarize


def packet(fiber, ra, rest=(3934.77,), z=.04, direction="absorption"):
    obs = f"sdss-6497-56329-{fiber:04d}"
    features = [{"classification": "unidentified", "direction": direction, "score": 8.0,
                 "rest_wavelength": wave, "wavelength": wave * (1 + z)} for wave in rest]
    return {"schema_version": 1, "state": "completed", "claim_status": "screening_only",
            "software_sha256": "a" * 64,
            "observation": {"observation_id": obs, "source_sha256": f"{fiber:064x}", "fiber": fiber,
                            "ra_deg": ra, "dec_deg": 0., "plate": 6497, "mjd": 56329,
                            "wavelength_unit": "Angstrom", "wavelength_frame": "observed vacuum"},
            "stages": {"detection": {"schema_version": 1, "score_threshold": 6., "method": "method-1",
                                     "catalog": "lines-1", "features": features,
                                     "source_context": {"kind": "observed", "observation_id": obs}}}}


class RecurrenceTests(unittest.TestCase):
    def test_measured_k_case_and_fixed_observer_comparison(self):
        # Frozen centroids and redshifts from the pipeline's three published coadds.
        sources = [packet(88, 214.51434, (3933.8294823692618,), .0422869473695755),
                   packet(96, 214.37598, (3934.71044975572,), .03487289324402809),
                   packet(100, 214.44446, (3935.4351489630185,), .07697335630655289)]
        result = summarize(sources)
        group, = result["rest_frame_groups"]
        self.assertEqual(group["target_count"], 3)
        self.assertEqual(group["visit_count"], 1)
        self.assertAlmostEqual(group["observed_span_angstrom"], 166.4336, places=3)
        self.assertEqual(result["observer_frame_control_groups"], [])
        self.assertEqual(result, summarize(reversed(sources)))

    def test_compatible_reviews_duplicate_sources_and_incomplete_packets(self):
        p = packet(1, 0.)
        reviewed = copy.deepcopy(p)
        reviewed["stages"]["reviewer"] = {"summary": "different model"}
        pending = {**p, "state": "queued", "stages": {}}
        result = summarize([p, reviewed, pending])
        self.assertEqual(result["unique_observations"], 1)
        self.assertEqual(result["duplicate_completed_packets"], 1)
        self.assertEqual(result["skipped_incomplete"], 1)
        self.assertEqual(result["rest_frame_groups"], [])
        conflicted = copy.deepcopy(p)
        conflicted["stages"]["detection"]["features"][0]["score"] = 10.
        with self.assertRaisesRegex(ValueError, "Conflicting evidence"):
            summarize([p, conflicted])

    def test_incompatible_settings_rejected(self):
        for location, field in (("packet", "software_sha256"), ("scan", "catalog"), ("scan", "method"),
                                ("scan", "schema_version"), ("scan", "score_threshold")):
            a, b = packet(1, .1), packet(2, 3.)
            target = b if location == "packet" else b["stages"]["detection"]
            target[field] = "changed"
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "Incompatible scans"):
                summarize([a, b])

    def test_ra_wrap_repeat_visits_and_one_vote_per_target(self):
        east = packet(1, 359.9999, (5000., 5000.3))
        west = packet(2, .0001, (5000.2,))
        west["observation"]["mjd"] += 1
        far = packet(3, 3.0, (5000.1,))
        self.assertEqual(summarize([east, west])["rest_frame_groups"], [])
        result = summarize([west, far, east])
        self.assertEqual(result["distinct_targets"], 2)
        group, = result["rest_frame_groups"]
        self.assertEqual(len(group["members"]), 2)
        self.assertEqual(result, summarize([east, far, west]))

    def test_pairs_cannot_chain_into_an_out_of_tolerance_triple(self):
        sources = [packet(i, float(i), (wave,)) for i, wave in enumerate((5000., 5001.8, 5003.6), 1)]
        groups = summarize(sources)["rest_frame_groups"]
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["target_count"], 2)
        self.assertLessEqual(groups[0]["frame_span_kms"], 150.)

    def test_sign_classification_and_observer_frame_control(self):
        a = packet(1, .1, (5000,), z=.1)
        b = packet(2, 2., (4400,), z=.25)
        c = packet(3, 3., (5000,), direction="emission")
        d = packet(4, 4., (5000,))
        d["stages"]["detection"]["features"][0]["classification"] = "catalog_match"
        result = summarize([a, b, c, d])
        self.assertEqual(result["rest_frame_groups"], [])
        control, = result["observer_frame_control_groups"]
        self.assertEqual(control["observed_span_angstrom"], 0.)
        self.assertEqual(result["unidentified_features"], 3)

    def test_resolve_only_campaign_reports_or_explicit_files(self):
        with tempfile.TemporaryDirectory() as name:
            campaign = Path(name)
            reportdir = campaign / "pages" / "first" / "reports"
            reportdir.mkdir(parents=True)
            evidence = reportdir / "example.json"
            evidence.write_text(json.dumps(packet(1, 0.)))
            (campaign / "state.json").write_text("{}")
            self.assertEqual(packet_files([campaign, reportdir, evidence]), [evidence.resolve()])


if __name__ == "__main__":
    unittest.main()
