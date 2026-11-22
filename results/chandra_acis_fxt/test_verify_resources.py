"""The runner-memory estimate can never skip execution of the repository verifier."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import verify_resources as budget


class MemoryBudget(unittest.TestCase):
    tree={"by_tracked_path":{"all":{"bytes":4*1024**3}}}
    low={"proc_meminfo":{"MemAvailable":"8000000 kB","SwapFree":"0 kB"},"disk_free_bytes":100*1024**3}

    def test_capacity_and_scope(self):
        generous={**self.low,"proc_meminfo":{"MemAvailable":"30000000 kB","SwapFree":"2000000 kB"}}
        self.assertEqual(budget.reserve_for_tree(self.tree,generous)["status"],"sufficient")
        with patch.dict(os.environ,{},clear=True):
            self.assertEqual(budget.reserve_for_tree(self.tree,self.low)["status"],"outside_research_runner")

    def test_unavailable_reservation_leaves_no_file(self):
        with tempfile.TemporaryDirectory() as folder:
            env={"GITHUB_ACTIONS":"true","GITHUB_REF":"refs/heads/research/chandra-acis-fast-transients",
                 "GITHUB_WORKFLOW_REF":"dovvnloading/Open-Space/.github/workflows/chandra-acis-fxt.yml@refs/heads/research/chandra-acis-fast-transients",
                 "RUNNER_TEMP":folder}
            fail=subprocess.CompletedProcess([],1,"unavailable")
            with patch.dict(os.environ,env), patch.object(budget.subprocess,"run",return_value=fail) as command:
                self.assertEqual(budget.reserve_for_tree(self.tree,self.low)["status"],"unavailable")
                command.assert_called_once()
                self.assertEqual(list(Path(folder).iterdir()),[])
                self.assertEqual(budget.reserve_for_tree(self.tree,{**self.low,"disk_free_bytes":1})["status"],"insufficient_free_disk")
                command.assert_called_once()
