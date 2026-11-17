"""Offline regressions for target-bound, frozen NEXUS stamp caches."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

MODULE = Path(__file__).resolve().parents[1] / 'results/jwst-repeat-imaging-transients/high-z-hosts/acquire_nexus.py'
spec = importlib.util.spec_from_file_location('nexus_cache_replay', MODULE)
acquire = importlib.util.module_from_spec(spec)
spec.loader.exec_module(acquire)


class FrozenStampCache(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cache = Path(self.temp.name)
        self.targets = pd.DataFrame([dict(id=1, ra=1., dec=2.)])
        self.npz = self.cache / 'ep01_data.npz'
        self.meta = self.cache / 'ep01_data.json'
        np.savez_compressed(self.npz, stamps=np.zeros((1, 64, 64)), px=[0.], py=[0.])
        frozen = next(row for row in json.loads((acquire.HERE / 'f200w/mosaic_inputs.json').read_text())
                      if row['epoch'] == 1 and row['kind'] == 'data')
        self.record = dict(frozen, cache_identity=dict(
            targets_sha256=hashlib.sha256(self.targets.to_csv(index=False).encode()).hexdigest(),
            radius=acquire.R, seed=acquire.SEED, binning=2,
            source_sha256=frozen['sha256'], source_bytes=frozen['bytes']),
            stamp_sha256=hashlib.sha256(self.npz.read_bytes()).hexdigest())
        self.write_record(self.record)

    def write_record(self, record):
        self.meta.write_text(json.dumps(record))

    def extract(self, targets=None):
        with patch.object(acquire, 'urlopen', side_effect=AssertionError('unexpected network')):
            acquire.extract((1, 'data', self.targets if targets is None else targets, self.cache, 'f200w'))

    def test_matching_cache_is_reused_without_changes(self):
        before = (self.npz.read_bytes(), self.meta.read_bytes())
        self.extract()
        self.assertEqual((self.npz.read_bytes(), self.meta.read_bytes()), before)

    def test_changed_target_list_is_rejected_without_writes(self):
        before = (self.npz.read_bytes(), self.meta.read_bytes())
        changed = self.targets.copy()
        changed.loc[0, 'id'] = 2
        with self.assertRaises(ValueError):
            self.extract(changed)
        self.assertEqual((self.npz.read_bytes(), self.meta.read_bytes()), before)

    def test_header_and_source_metadata_mutants_are_rejected(self):
        for field in ('header', 'sha256', 'bytes', 'filter', 'binning'):
            with self.subTest(field=field):
                changed = copy.deepcopy(self.record)
                if field == 'header':
                    changed['header']['PIXAR_SR'] = 999
                else:
                    changed[field] = 'changed'
                self.write_record(changed)
                before = (self.npz.read_bytes(), self.meta.read_bytes())
                with self.assertRaises(ValueError):
                    self.extract()
                self.assertEqual((self.npz.read_bytes(), self.meta.read_bytes()), before)

    def test_changed_stamp_bytes_and_partial_cache_are_rejected(self):
        self.npz.write_bytes(b'changed payload')
        with self.assertRaises(ValueError):
            self.extract()
        self.assertEqual(self.npz.read_bytes(), b'changed payload')
        self.meta.unlink()
        with self.assertRaisesRegex(ValueError, 'Partial stamp cache'):
            self.extract()
        self.assertFalse(self.meta.exists())


if __name__ == '__main__':
    unittest.main()
