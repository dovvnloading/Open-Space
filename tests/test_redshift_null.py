import copy
import unittest

from open_space.redshift_null import _assignments, experiment
from open_space.spectra import CATALOG_VERSION, _identify


def packet(i, z, rest, ra=None):
    obs_id = f"sample-{i}"
    obs = {"observation_id": obs_id, "source_sha256": f"{i:064x}", "csv_sha256": f"{i+10:064x}",
           "plate": 5, "mjd": 55555, "ra_deg": i if ra is None else ra, "dec_deg": 0,
           "wavelength_frame": "observed vacuum", "wavelength_unit": "Angstrom", "redshift": z}
    feats = [_identify({"direction": "absorption", "score": 8.0, "pixel_spacing_angstrom": 1.,
                       "wavelength": wavelength * (1 + z)}, z) for wavelength in rest]
    return {"schema_version": 1, "claim_status": "screening_only", "state": "completed",
            "software_sha256": "a" * 64, "observation": obs, "stages": {"detection": {
                "source_context": {"kind": "observed", "observation_id": obs_id},
                "features": feats, "method": "unit-test", "redshift": z, "score_threshold": 6.,
                "schema_version": 1, "catalog": CATALOG_VERSION}}}


class LabelShuffleTests(unittest.TestCase):
    def test_hand_enumeration_and_catalog_membership(self):
        # Swapping redshifts splits the pair of unknown lines. It also releases
        # both known H-beta detections when source matching is recomputed.
        packets = [packet(1, .03, [4500., 4862.683]), packet(2, .08, [4500., 4862.683])]
        frozen = copy.deepcopy(packets)
        r = experiment(packets, draws=2)
        self.assertEqual(r, experiment(reversed(packets), draws=2))
        self.assertEqual(packets, frozen)
        self.assertEqual(r["unrestricted"]["mode"], "exact")
        d = r["unrestricted"]["distributions"]
        self.assertEqual(d["reapply_catalog"]["features_after_catalog_screen"]["histogram"], {"2": 1, "4": 1})
        self.assertEqual(d["freeze_original_shortlist"]["features_after_catalog_screen"]["histogram"], {"2": 2})
        maximum = d["reapply_catalog"]["largest_group"]
        self.assertEqual(maximum["histogram"], {"0": 1, "2": 1})
        self.assertEqual(maximum["upper_tail_fraction"], .5)
        self.assertEqual(r["observed"]["largest_group"], 2)

    def test_whole_spectra_and_zero_feature_targets(self):
        r = experiment([packet(1, .03, [4500.]), packet(2, .05, [4500.]), packet(3, .08, [])], draws=10)
        self.assertEqual(r["unrestricted"]["possible_index_assignments"], 6)
        self.assertEqual(r["unrestricted"]["assignments"], 6)
        self.assertEqual(r["observed"]["largest_group"], 2)
        self.assertEqual(len(r["inputs"]), 3)

    def test_monte_carlo_tail_and_neighbor_block_membership(self):
        sample = [packet(i, .02 * i, [4567.]) for i in range(1, 7)]
        r = experiment(sample, draws=8, seed=41)
        e = r["unrestricted"]["distributions"]["reapply_catalog"]["largest_group"]
        self.assertEqual(r["unrestricted"]["mode"], "monte_carlo_with_replacement")
        self.assertEqual(e["upper_tail_fraction"], (1 + e["draws_at_least_observed"]) / 9)
        self.assertEqual(r["neighboring_redshift_triples"]["possible_index_assignments"], 36)
        _, it = _assignments([[0, 2, 4], [1, 3, 5]], 36, 0)
        all_assignments = {tuple(p) for p in it}
        self.assertEqual(len(all_assignments), 36)
        for p in all_assignments:
            self.assertEqual(set(p[::2]), {0, 2, 4})
            self.assertEqual(set(p[1::2]), {1, 3, 5})

    def test_identity_failure_and_invalid_input(self):
        a, b = packet(1, .03, [4500.]), packet(2, .08, [4500.])
        for bad in ([a, a], [a, packet(2, .08, [4500.], ra=1.)],
                    [a, {**b, "state": "running"}], [a, packet(2, .03, [4500.])]):
            with self.assertRaises(ValueError):
                experiment(bad, draws=2)
        corrupt = copy.deepcopy(b)
        corrupt["stages"]["detection"]["features"][0]["classification"] = "catalog_match"
        with self.assertRaisesRegex(ValueError, "Identity"):
            experiment([a, corrupt], draws=2)
        corrupt = copy.deepcopy(b)
        corrupt["observation"]["mjd"] += 1
        with self.assertRaisesRegex(ValueError, "one plate/MJD"):
            experiment([a, corrupt], draws=2)


if __name__ == "__main__":
    unittest.main()
