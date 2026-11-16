"""The repository must pass its own verification.

This is the gate: a pull request that breaks a checksum, a link, the catalog
counts or a declared claim fails here, in the ordinary test run.
"""

import shutil
import unittest
from pathlib import Path

from open_space.verify import runner
from open_space.verify.model import FAIL
from open_space.verify.tree import RepoTree, is_git_repo

ROOT = Path(__file__).resolve().parents[1]


class RepositoryVerifiesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("git") and is_git_repo(ROOT):
            tree = RepoTree.from_worktree(ROOT)          # the tree as it would be committed
        else:
            tree = RepoTree.from_directory(ROOT)
        cls.ledger = runner.run_checks(tree)
        cls.summary = cls.ledger["summary"]

    def test_no_machine_checkable_statement_is_false(self):
        failures = [f"{finding.rule} at {finding.where()}: {finding.message}"
                    for finding in self.ledger["findings"] if finding.severity == FAIL]
        self.assertEqual(failures, [], "\n" + "\n".join(failures))

    def test_catalog_counts_and_links_hold(self):
        counts = self.summary["catalog_counts"]
        self.assertEqual(counts["stated_total"], counts["entries"])
        for label, stated in counts["stated"].items():
            self.assertEqual(stated, counts["labelled_entries"][label], label)
        self.assertEqual(self.summary["links"]["broken"], 0)

    def test_recorded_checksums_match_committed_bytes(self):
        checksums = self.summary["checksums"]
        self.assertEqual(checksums["mismatched"], 0)
        self.assertGreater(checksums["verified_against_committed_bytes"], 1000)

    def test_every_declared_claim_recomputes(self):
        claims = self.summary["claims"]
        self.assertGreaterEqual(claims["claims"], 50)
        self.assertEqual(claims["reproduced"], claims["claims"])

    def test_catalog_numbers_appear_in_the_results_they_cite(self):
        self.assertEqual(self.summary["numbers"]["catalog"]["unmatched"], 0)


if __name__ == "__main__":
    unittest.main()
