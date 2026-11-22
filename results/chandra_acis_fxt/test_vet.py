"""Offline fixtures for catalog classification and exact-response replay."""
import hashlib
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from xml.sax.saxutils import escape

import vet


TRIGGER = {"ra_deg": 192.0, "dec_deg": 27.0, "epoch_year": 2024.5,
           "position_error_arcsec": 1.5, "obsid": 19000}
CUDIDI_HEADER = b"CANDIDATE,OBSID,RA,DEC\n"


def votable(columns, row=(), repetitions=0, trailing_status=""):
    fields = "".join(f'<FIELD name="{col}" datatype="char" arraysize="*"/>' for col in columns)
    cells = "".join(f"<TD>{escape(str(v))}</TD>" for v in row)
    rows = f"<TR>{cells}</TR>" * repetitions
    trailer = f'<INFO name="QUERY_STATUS" value="{trailing_status}"/>' if trailing_status else ""
    return (f'<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">'
            f'<RESOURCE type="results"><INFO name="QUERY_STATUS" value="OK"/><TABLE>{fields}'
            f'<DATA><TABLEDATA>{rows}</TABLEDATA></DATA></TABLE>{trailer}</RESOURCE></VOTABLE>').encode()


class TestVet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bodies = {
            "gaia": votable("source_id ra dec ref_epoch pmra pmdec pmra_error pmdec_error "
                            "ra_error dec_error parallax parallax_error ruwe phot_g_mean_mag".split()),
            "simbad": votable("main_id ra dec otype rvz_redshift path is_candidate coo_err_maj".split()),
            "ned": votable("prefname ra dec pretype z uncmaja".split()),
            # Published whole tables have fixed sizes. All fixture positions are far away.
            "quirola2023": votable("XRT RAdeg DEdeg ObsID ePos".split(),
                                   ["distant-fixture", 0, 0, 10, 1], repetitions=8),
            "dillmann2025": votable("CSCName RACdeg DECdeg CSCobsId TransType".split(),
                                    ["distant-fixture", 0, 0, 11, "F"], repetitions=3559),
            "cudidi2026": CUDIDI_HEADER + b"far-fixture,12,0.0,0.0\n",
        }

    def screen(self, trigger=None, replacements=None):
        bodies = {**self.bodies, **(replacements or {})}

        def provide(label, url, root, replay, records, params=None):
            self.assertTrue(replay)
            body = bodies[label]
            records[label] = {"url": url, "sha256": hashlib.sha256(body).hexdigest()}
            return body

        # Substitute the checksum published by this synthetic fixture. The normal
        # production decoder, VOTable validation and classification paths still run.
        csv_md5 = hashlib.md5(bodies["cudidi2026"], usedforsecurity=False).hexdigest()
        with tempfile.TemporaryDirectory() as work, patch.object(vet, "read_input", side_effect=provide), \
                patch.object(vet, "CUDIDI_MD5", csv_md5), patch.object(vet, "urlopen") as network:
            result = vet.vet_trigger(trigger or TRIGGER, work, replay=True)
            network.assert_not_called()
            return result

    def test_empty_success_and_truncated_query_have_different_dispositions(self):
        empty = self.screen()
        self.assertEqual(empty["status"], "no_counterpart")
        self.assertEqual(empty["matches"], [])
        self.assertFalse(empty["errors"])
        self.assertEqual(empty["queries"]["ned"]["row_count"], 0)

        truncated = votable("prefname ra dec pretype z uncmaja".split(), trailing_status="OVERFLOW")
        incomplete = self.screen(replacements={"ned": truncated})
        self.assertEqual(incomplete["status"], "catalog_incomplete")
        self.assertIn("ned", incomplete["errors"])
        self.assertEqual(incomplete["queries"]["ned"]["sha256"], hashlib.sha256(truncated).hexdigest())

    def test_published_position_in_same_and_different_visit(self):
        fixture = {"cudidi2026": CUDIDI_HEADER + b"XT fixture,19000,192.0,27.0\n"}
        published = self.screen(replacements=fixture)
        self.assertEqual(published["status"], "known_transient")
        self.assertEqual(published["matches"][0]["name"], "XT fixture")
        other_visit = self.screen(trigger={**TRIGGER, "obsid": 19001}, replacements=fixture)
        self.assertEqual(other_visit["status"], "ambiguous")
        self.assertEqual(other_visit["matches"][0]["kind"], "published_variability")

    def test_mirror_preserves_failed_and_successful_receipts(self):
        unavailable = b"<html>Service unavailable</html>"
        result = self.screen(replacements={"gaia": unavailable, "gaia_aip": self.bodies["gaia"]})
        self.assertEqual(result["status"], "no_counterpart")
        self.assertIn("query_error", result["queries"]["gaia"])
        self.assertEqual(result["queries"]["gaia"]["sha256"], hashlib.sha256(unavailable).hexdigest())
        self.assertEqual(result["queries"]["gaia_aip"]["row_count"], 0)

    def test_warm_cache_and_replay_verify_bytes_without_network(self):
        original = b"fixture A"
        url = "https://catalog.example.invalid/fixture"

        class Response(io.BytesIO):
            status = 200

            def geturl(self):
                return url

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = {}
            workflow = {"GITHUB_ACTIONS": "true", "GITHUB_EVENT_NAME": "push", "GITHUB_REF": vet.CAMPAIGN_REF}
            with patch.dict(os.environ, workflow), patch.object(vet, "urlopen", return_value=Response(original)) as network:
                self.assertEqual(vet.read_input("fixture", url, root, False, records), original)
                self.assertEqual(network.call_count, 1)
            with patch.object(vet, "urlopen", side_effect=AssertionError("cache replay contacted network")) as network:
                self.assertEqual(vet.read_input("fixture", url, root, False, {}), original)
                self.assertEqual(vet.read_input("fixture", url, root, True, {}), original)
                (root / records["fixture"]["body_file"]).write_bytes(b"fixture B")
                with self.assertRaisesRegex(ValueError, "altered catalog response"):
                    vet.read_input("fixture", url, root, True, {})
                network.assert_not_called()

    def test_push_on_another_branch_cannot_fetch(self):
        env = {"GITHUB_ACTIONS": "true", "GITHUB_EVENT_NAME": "push", "GITHUB_REF": "refs/heads/main"}
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, env), patch.object(vet, "urlopen") as network:
            with self.assertRaisesRegex(RuntimeError, "research branch"):
                vet.read_input("fixture", "https://catalog.example.invalid/fixture", Path(tmp), False, {})
            network.assert_not_called()
