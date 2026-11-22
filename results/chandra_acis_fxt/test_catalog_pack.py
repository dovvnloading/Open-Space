"""Synthetic regression checks for the optional final evidence packer."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

import catalog_pack as cp


class CatalogPackTest(unittest.TestCase):
    def setUp(self):
        area = tempfile.TemporaryDirectory()
        self.addCleanup(area.cleanup)
        self.area = Path(area.name)
        self.root = self.area / "original"
        self.ns = "run-1"
        self.folder = self.root / "catalogs" / self.ns
        self.folder.mkdir(parents=True)
        self.files = {
            "visible.bin": b"retained for written review",
            "empty.bin": b"",
            "small.bin": bytes(range(251)) * 4,
            "big.bin": bytes(range(251)) * 48,
        }
        for name, body in self.files.items():
            (self.folder / name).write_bytes(body)
            entry = {"body_file": name, "size_bytes": len(body), "sha256": cp.digest(body)}
            memo = {"url": "https://example.invalid/query", "attempts": [entry], **entry}
            (self.folder / (name + ".json")).write_text(json.dumps(memo))
        self.keep = {self.logical("visible.bin")}
        for field, value in (("BLOCK", 4096), ("PACK", 8192), ("MAX_MEMBERS", 3)):
            mock = patch.object(cp, field, value)
            mock.start()
            self.addCleanup(mock.stop)

    def logical(self, name):
        return (self.folder / name).relative_to(self.root).as_posix()

    def store(self):
        return cp.Store(self.root / "catalog_store" / self.ns, self.ns)

    def test_round_trip_is_deterministic_and_preserves_linked_and_empty_sources(self):
        proposed = cp.pack(self.root, self.ns, self.keep, True, False)
        self.assertFalse((self.root / "catalog_store").exists())
        self.assertEqual(set(self.files), {p.name for p in self.folder.glob("*.bin")})
        published = cp.pack(self.root, self.ns, self.keep, False, False)
        self.assertEqual(proposed, published)
        before = (self.store().dir / "index.json").read_bytes()
        self.assertEqual(set(self.files), {p.name for p in self.folder.glob("*.bin")})
        cp.pack(self.root, self.ns, self.keep, False, True)
        self.assertEqual({"visible.bin"}, {p.name for p in self.folder.glob("*.bin")})
        self.assertEqual(4, len(list(self.folder.glob("*.json"))))
        cp.pack(self.root, self.ns, self.keep, False, True)
        db = self.store()
        self.assertEqual(before, (db.dir / "index.json").read_bytes())
        self.assertGreater(len(db.entries[self.logical("big.bin")]), 1)
        self.assertEqual(1, len(db.entries[self.logical("empty.bin")]))
        for pack in db.index["packs"]:
            self.assertEqual(pack["sha256"], cp.sha(db.dir / pack["file"]))
            with ZipFile(db.dir / pack["file"]) as z:
                for member in z.infolist():
                    self.assertLessEqual(member.file_size, cp.BLOCK)
        cp.verify_sources(db, self.root)
        view = self.area / "view"
        cp.extract(self.root, self.ns, view)
        for name, body in self.files.items():
            self.assertEqual(body, (view / self.logical(name)).read_bytes())
        (view / self.logical("visible.bin")).write_bytes(b"changed privately")
        self.assertEqual(self.files["visible.bin"], (self.folder / "visible.bin").read_bytes())
        cp.pack(self.root, self.ns, self.keep | {self.logical("big.bin")}, False, True)
        self.assertEqual(self.files["big.bin"], (self.folder / "big.bin").read_bytes())
        cp.verify_sources(self.store(), self.root)

    def test_corrupt_pack_is_rejected_and_raw_replies_stay_available(self):
        cp.pack(self.root, self.ns, self.keep, False, False)
        pack = next(self.store().dir.glob("*.zip"))
        with pack.open("ab") as f:
            f.write(b"corruption")
        with self.assertRaisesRegex(ValueError, "content address"):
            cp.pack(self.root, self.ns, self.keep, False, True)
        self.assertEqual(set(self.files), {p.name for p in self.folder.glob("*.bin")})

    def test_selective_read_opens_only_the_bundles_spanning_the_reply(self):
        cp.pack(self.root, self.ns, self.keep, False, True)
        store = self.store()
        name = "small.bin"
        relevant = {pack for pack, _ in store.entries[self.logical(name)]}
        unrelated = next(p for p in store.dir.glob("*.zip") if p.name not in relevant)
        unrelated.unlink()
        body = self.folder / name
        receipt = json.loads((self.folder / (name + ".json")).read_text())
        private = self.area / "private" / name
        cp.recorded_reply(self.root, body, receipt, private)
        self.assertEqual(private.read_bytes(), self.files[name])
        self.assertFalse(body.exists())
        with self.assertRaises(OSError):
            self.store()  # Final whole-store verification still reports the missing unrelated pack.

    def test_changed_or_missing_memo_is_refused(self):
        memo = self.folder / "small.bin.json"
        original = memo.read_bytes()
        data = json.loads(original)
        data["attempts"][0]["sha256"] = "0" * 64
        memo.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "memo/body"):
            cp.pack(self.root, self.ns, self.keep, False, True)
        self.assertTrue((self.folder / "small.bin").exists())
        memo.write_bytes(original)
        cp.pack(self.root, self.ns, self.keep, False, False)
        memo.unlink()
        with self.assertRaisesRegex(ValueError, "removed or altered"):
            cp.pack(self.root, self.ns, self.keep, False, True)
        self.assertTrue((self.folder / "small.bin").exists())


if __name__ == "__main__":
    unittest.main()
