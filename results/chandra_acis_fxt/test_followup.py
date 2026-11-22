"""Network-isolated regressions for the selective repair and its evidence history."""
from copy import deepcopy
import hashlib
import io
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import campaign
import followup as fix
import test_vet as fixtures

TRIGGER = fixtures.TRIGGER
votable = fixtures.votable


def ledger(oid, vetting=None):
    h = hashlib.sha256(f"synthetic event {oid}".encode()).hexdigest()
    product = {"filename": f"acisf{oid:05d}N001_evt2.fits.gz", "sha256": h, "bytes": 1,
               "obsid": oid, "url": "https://cxc.cfa.harvard.edu/fixture"}
    seed = int(hashlib.sha256(("chandra-acis-fxt-v1:" + h).encode()).hexdigest()[:16], 16) >> 1
    stream = {"seed": seed, "obsid": oid, "triggers": [],
              "controls": [{"replica": n, "triggers": []} for n in range(24)],
              "injections": [{"trial": n} for n in range(12)]}
    if vetting is not None:
        for key, events in (("real/0", stream["triggers"]), ("control-0/0", stream["controls"][0]["triggers"])):
            events.append({**TRIGGER, "obsid": oid, "trigger_key": key, "vet": deepcopy(vetting),
                           "instrument_pass": False, "catalog_cache": "catalogs/primary"})
    return {"obsid": oid, "status": "complete", "catalog_cache": "catalogs/primary", "errors": [],
            "acquired": [product], "products": [{"input": product, "complete": True, "search": stream}]}


class Response(io.BytesIO):
    status = 200

    def geturl(self):
        return "https://fixture.invalid/recorded"


class TestFollowup(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.work = Path(tmp.name)
        self.evidence, self.bundle = self.work / "evidence", self.work / "bundle"
        (self.evidence / "observations").mkdir(parents=True)
        # The reducer may install reviewed source, always inside this isolated copy.
        source = self.work / "implementation"
        source.mkdir()
        for name in (*campaign.identity(), "replay_untyped_ned.py"):
            shutil.copyfile(campaign.HERE / name, source / name)
        self.enterContext(patch.object(campaign, "HERE", source))
        self.enterContext(patch.dict(os.environ, {"GITHUB_ACTIONS": "true",
           "GITHUB_REF": "refs/heads/research/chandra-acis-fast-transients", "GITHUB_EVENT_NAME": "workflow_dispatch",
           "GITHUB_RUN_ID": "12345", "GITHUB_RUN_ATTEMPT": "1", "GITHUB_SHA": "fixture",
           "GITHUB_OUTPUT": str(self.work / "output")}))
        original = sys.modules.get("vet")
        self.addCleanup(lambda: sys.modules.__setitem__("vet", original))
        code = self.work / "source-pair"
        fix.source_pair(code)
        self.old = fix.load_vet(code / "vet.before.py")
        self.after_path = code / "vet.after.py"
        # No fixture can make a real remote request, including the original vet.
        self.enterContext(patch.object(self.old, "urlopen", side_effect=AssertionError("unexpected external request")))

    def historic(self, primary_ned, mirror_ned):
        fixtures.TestVet.setUpClass()
        bodies = {**fixtures.TestVet.bodies, "ned": primary_ned, "ned_vo": mirror_ned}
        self.fixture_md5 = hashlib.md5(bodies["cudidi2026"], usedforsecurity=False).hexdigest()
        reader = self.old.read_input
        def read(label, url, root, replay, records, params=None):
            with patch.object(self.old, "urlopen", return_value=Response(bodies[label])):
                return reader(label, url, root, replay, records, params)
        folder = self.evidence / "catalogs/primary"
        with patch.object(self.old, "CUDIDI_MD5", self.fixture_md5), patch.object(self.old, "read_input", read):
            return self.old.vet_trigger(TRIGGER, folder)

    def corrected(self):
        corrected = fix.load_vet(self.after_path)
        self.enterContext(patch.object(corrected, "CUDIDI_MD5", self.fixture_md5))
        self.enterContext(patch.object(corrected, "urlopen", side_effect=AssertionError("offline replay contacted network")))
        return corrected

    @staticmethod
    def ned(empty=False):
        return votable("prefname ra dec preftype z uncmaja".split(), ["untyped", 192, 27, "", 0, 1],
                       repetitions=0 if empty else 1)

    def check_blank(self, mirror_empty):
        v = self.historic(self.ned(), self.ned(empty=mirror_empty))
        self.assertTrue(fix.affected(v))
        if mirror_empty:
            self.assertEqual(v["status"], "no_counterpart")
            self.assertEqual(v["queries"]["ned"]["query_error"], fix.UNKNOWN)
        else:
            self.assertEqual(v["errors"]["ned"], fix.UNKNOWN)
        row = ledger(19000, v)
        p = row["products"][0]
        before = deepcopy(p)
        prior = {f.name: fix.sha(f) for f in (self.evidence / "catalogs/primary").iterdir()}
        corrected = self.corrected()
        out, changes = self.work / "worker", {}
        for stream, i, event in fix.triggers(p):
            fix.prepare_event(row, p, stream, i, event, self.evidence, out, corrected, {}, changes)
            self.assertEqual(event["vet"]["status"], "ambiguous")
            self.assertEqual([(m["catalog"], m["kind"]) for m in event["vet"]["matches"]],
                             [("ned", "unclassified")])
            self.assertEqual(event["vet"]["queries"]["ned"]["sha256"], v["queries"]["ned"]["sha256"])
            self.assertNotIn("ned_vo", event["vet"]["queries"])
            self.assertTrue((out / event["catalog_cache"] / v["queries"]["ned"]["body_file"]).is_file())
        self.assertEqual(len(changes), 2)
        self.assertEqual(fix.raw_scan(p), fix.raw_scan(before))
        self.assertEqual(prior, {f.name: fix.sha(f) for f in (self.evidence / "catalogs/primary").iterdir()})

    def test_all_blank_real_and_control_replay_exact_response(self):
        self.check_blank(mirror_empty=False)

    def test_previously_successful_mirror_does_not_hide_blank_primary(self):
        self.check_blank(mirror_empty=True)

    def test_unrelated_primary_failure_with_successful_mirror_is_unchanged(self):
        v = self.historic(b"<html>Unavailable</html>", self.ned(empty=True))
        self.assertTrue(fix.certify(v))
        self.assertFalse(fix.affected(v))
        self.assertEqual(v["queries"]["ned_vo"]["sha256"], hashlib.sha256(self.ned(empty=True)).hexdigest())
        row = ledger(19000, v)
        self.assertEqual(fix.inspect(19000, row, self.evidence)[0], [])
        (self.evidence / "catalogs/primary" / v["queries"]["ned_vo"]["body_file"]).unlink()
        self.assertIn("missing_or_altered_original_catalog_inputs", fix.inspect(19000, row, self.evidence)[0])
        with self.assertRaisesRegex(ValueError, "missing or altered"):
            fix.seed_cache(self.evidence, row, row["products"][0]["search"]["triggers"][0],
                           self.work / "retry", live=True)

    def test_missing_cache_memo_can_be_rebuilt_from_original_raw_bytes(self):
        v = self.historic(b"<html></html>", self.ned(empty=True))
        q = v["queries"]["ned_vo"]
        key = hashlib.sha256(q["url"].encode()).hexdigest()[:24]
        (self.evidence / "catalogs/primary" / f"ned_vo-{key}.json").unlink()
        row = ledger(19000, v)
        self.assertIn("catalog_memo_requires_reconstruction", fix.inspect(19000, row, self.evidence)[0])
        p, changed = row["products"][0], {}
        stream, i, t = next(fix.triggers(p))
        fix.prepare_event(row, p, stream, i, t, self.evidence, self.work / "worker", self.corrected(), {}, changed)
        self.assertEqual(t["vet"], v)
        self.assertTrue(fix.memo_matches(self.work / "worker" / t["catalog_cache"], v))

    def test_packed_successes_remain_equivalent_and_detect_corruption(self):
        v = self.historic(b"<html>Unavailable</html>", self.ned(empty=True))
        row = ledger(19000, v)
        expected = fix.inspect(19000, row, self.evidence)
        self.assertEqual(expected[0], [])
        unrelated = self.evidence / "catalog_store" / "unrelated"
        unrelated.mkdir(parents=True)
        (unrelated / "index.json").write_text("This intentionally invalid unrelated cache must stay unopened.")
        fix.catalog_pack.pack(self.evidence, "primary", set(), False, True)
        self.assertFalse(any((self.evidence / "catalogs/primary").glob("*.bin")))
        self.assertEqual(expected, fix.inspect(19000, row, self.evidence))
        memo = next((self.evidence / "catalogs/primary").glob("ned_vo-*.json"))
        original = memo.read_bytes()
        memo.write_bytes(original + b"\n")
        self.assertIn("catalog_memo_requires_reconstruction", fix.inspect(19000, row, self.evidence)[0])
        memo.write_bytes(original)
        # A prior verification in this same process must not conceal changed ZIP bytes.
        part = next((self.evidence / "catalog_store/primary").glob("*.zip"))
        with part.open("ab") as stream:
            stream.write(b"unexpected")
        self.assertIn("missing_or_altered_original_catalog_inputs", fix.inspect(19000, row, self.evidence)[0])

    def test_packed_and_loose_affected_real_and_control_replay_identically(self):
        v = self.historic(self.ned(), self.ned())
        base = ledger(19000, v)
        primary_cache = self.evidence / "catalogs/primary"
        (primary_cache / "unrelated.bin").write_bytes(b"Unrelated recorded body in this same producer")
        corrected = self.corrected()
        outputs = []
        for mode in ("loose", "packed"):
            if mode == "packed":
                fix.catalog_pack.pack(self.evidence, "primary", set(), False, True)
            row = deepcopy(base)
            product = row["products"][0]
            out = self.work / mode
            changes = {}
            for stream,i,event in fix.triggers(product):
                fix.prepare_event(row, product, stream, i, event, self.evidence, out, corrected, {}, changes)
                self.assertEqual(event["vet"]["status"], "ambiguous")
                stage = out / event["catalog_cache"]
                self.assertFalse((stage / "unrelated.bin").exists())
                fix.verify_original(event["vet"], stage, out)
            self.assertEqual(fix.inspect(19000, row, self.evidence, (out, fix.namespace()))[0], [])
            self.assertIn("missing_or_altered_original_catalog_inputs", fix.inspect(19000, row, self.evidence)[0])
            identities = {f.relative_to(out).as_posix():fix.sha(f) for f in out.rglob("*") if f.is_file()}
            outputs.append((row, changes, identities))
        self.assertEqual(*outputs)
        self.assertFalse(any(primary_cache.glob("*.bin")))
        self.assertTrue((self.evidence / "catalog_store/primary/index.json").is_file())

    def plan(self):
        fix.save(self.evidence / "inventory/inventory.json", [{"obsid": n, "event_count": n} for n in (1, 2, 3)])
        self.enterContext(patch.object(fix, "INVENTORY", fix.sha(self.evidence / "inventory/inventory.json")))
        freeze = {"implementation": {**campaign.identity(), "vet.py":fix.OLD}, "git_commit": "primary",
                  "n_controls":24, "n_injections_per_product":12, "inventory_sha256": fix.INVENTORY}
        fix.save(self.evidence / "inventory/campaign.json", freeze)
        for oid in (1, 2):
            row = ledger(oid)
            if oid == 2:
                row["status"] = "followup_incomplete"
            campaign.save_obs(fix.old_path(self.evidence, oid), row)
        fix.prepare(self.evidence, self.bundle, 2)
        return campaign.read(self.bundle / "plan.json")

    def test_absent_outputs_preserve_primary_and_resume_finalization(self):
        plan = self.plan()
        self.assertEqual({*sum(plan["groups"], [])}, {2, 3})
        primary = {oid: fix.old_path(self.evidence, oid).read_bytes() for oid in (1, 2)}
        empty = self.work / "empty"
        empty.mkdir()
        result = fix.finish(self.evidence, self.bundle, empty, True)
        self.assertEqual(result["remaining_obsids"], [2, 3])
        self.assertEqual(result["carried_complete_without_rescan"], 1)
        self.assertEqual(fix.old_path(self.evidence, 1).read_bytes(), primary[1])
        self.assertFalse(fix.old_path(self.evidence, 3).exists())
        history = self.evidence / fix.SCOPED / "12345-1"
        self.assertEqual(fix.old_path(history / "before", 2).read_bytes(), primary[2])
        self.assertEqual(fix.sha(history / "vet.before.py"), fix.OLD)
        self.assertEqual(fix.sha(history / "vet.after.py"), fix.NEW)
        self.assertEqual(campaign.read(self.evidence / "inventory/campaign.json")["implementation"], plan["after"])
        self.assertEqual(result, fix.finish(self.evidence, self.bundle, empty, False))
        later = self.work / "continuation"
        with patch.dict(os.environ, {"GITHUB_RUN_ID": "6789"}):
            fix.prepare(self.evidence, later, 2)
        continuation = campaign.read(later / "plan.json")
        self.assertEqual({*sum(continuation["groups"], [])}, {2, 3})
        self.assertEqual(continuation["before"], plan["after"])
        self.assertEqual(continuation["obsids"]["2"]["primary_sha256"], fix.sha(fix.old_path(self.evidence, 2)))

    def test_missing_controls_cannot_be_promoted(self):
        plan = self.plan()
        stage = self.work / "stage"
        row = ledger(3)
        row["products"][0]["search"]["controls"].pop()
        bad = fix.old_path(stage, 3)
        bad.parent.mkdir(parents=True)
        campaign.save_obs(bad, row)
        receipt = {"obsid":3, "outcome":"complete", "primary_sha256":None,
                   "revised_sha256":fix.sha(bad), "plan_sha256":fix.sha(self.bundle / "plan.json"),
                   "run_id":plan["run_id"], "run_attempt":plan["run_attempt"]}
        fix.save(stage / "repair_records/3.json", receipt)
        with self.assertRaisesRegex(ValueError, "claimed completion"):
            fix.finish(self.evidence, self.bundle, stage, True)
        self.assertFalse(fix.old_path(self.evidence, 3).exists())
        receipt["outcome"] = "unresolved"
        fix.save(stage / "repair_records/3.json", receipt)
        result = fix.finish(self.evidence, self.bundle, stage, True)
        self.assertIn(3, result["remaining_obsids"])
        revised = campaign.read_obs(fix.old_path(self.evidence, 3))
        self.assertEqual(revised["status"], "error")
        self.assertFalse(revised["products"][0]["complete"])

    def test_real_requires_photometry_and_all_controls_require_catalogs(self):
        v = self.historic(b"<html></html>", self.ned(empty=True))
        row = ledger(19000, v)
        p = row["products"][0]
        t = p["search"]["triggers"][0]
        t["instrument_pass"] = True
        self.assertIn("physical_characterization", fix.inspect(19000, row)[0])
        p["search"]["controls"][0]["triggers"][0]["vet"]["errors"] = {"gaia": "network error"}
        self.assertIn("catalog_requires_review", fix.inspect(19000, row)[0])
        t["flux"] = {"peak_flux_erg_cm2_s": 1e-14}
        self.assertIn("physical_characterization", fix.inspect(19000, row)[0])
        changed = deepcopy(row)
        changed["products"][0]["search"]["injections"][0]["recovered"] = True
        with self.assertRaisesRegex(ValueError, "injection outcomes changed"):
            fix.prior_preserved(row, changed)

    def test_event_retry_rejects_changed_bytes_and_keeps_mirror_receipt(self):
        original, changed = b"archived-fixture", b"revised-fixture"
        name = "acisf00002N001_evt2.fits.gz"
        p = {"filename": name, "obsid": 2, "url": fix.event_names(name, 2)[0],
             "bytes": len(original), "sha256": hashlib.sha256(original).hexdigest(), "fits_sha256": "f" * 64}
        sink = self.work / "retrieval"
        sink.mkdir()
        def download(url, path):
            raw = original if url == fix.event_names(name, 2)[1] else changed
            path.write_bytes(raw)
            return {"url": url, "sha256": hashlib.sha256(raw).hexdigest(), "bytes":len(raw)}
        with patch.object(fix.archive, "download", side_effect=download) as http, \
                patch.object(fix.archive, "inspect_product", return_value={"fits_sha256": "f" * 64}):
            trace = []
            product = fix.event_bytes(p, 2, sink, [], trace)
            self.assertEqual(http.call_count, 2)
            self.assertEqual(product["url"], p["url"])
            self.assertEqual(product["sha256"], p["sha256"])
            self.assertIn("differs", trace[0]["error"])
            self.assertEqual(trace[1]["response"]["url"], fix.event_names(name, 2)[1])
            self.assertEqual(trace[1]["fits_sha256"], p["fits_sha256"])
            self.assertEqual((sink / name).read_bytes(), original)
