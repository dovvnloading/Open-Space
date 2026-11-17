"""Regression tests for immutable catalog replay inputs."""
import csv
import hashlib
import importlib.util
from io import BytesIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / 'results/jwst-repeat-imaging-transients/repeat-catalogs/fetch_catalogs.py'
spec = importlib.util.spec_from_file_location('jwst_catalog_fetch', MODULE)
fetch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch)


class FrozenCatalogReplay(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root / 'cache'
        self.cache.mkdir()
        self.names = ['jw01208-o002_t001_nircam_clear-f115w', 'jw01208-o002_t001_nircam_clear-f150w']
        self.data = (' '.join(fetch.COLUMNS) + '\n1 10 20 3 1 0 20 3 3\n').encode()
        with (self.root / 'pairs.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['reference', 'search']); w.writeheader(); w.writerow(dict(zip(['reference', 'search'], self.names)))
        with (self.root / 'catalog_inputs.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['product', 'url', 'bytes', 'sha256', 'sources']); w.writeheader()
            for name in self.names:
                w.writerow(dict(product=name, url=fetch.product_url(name, 'cat.ecsv'), bytes=len(self.data), sha256=hashlib.sha256(self.data).hexdigest(), sources=1))
        self.receipt = (self.root / 'catalog_inputs.csv').read_bytes()

    def test_matching_download_and_reuse_rebuilds_untrusted_npz(self):
        with patch.object(fetch, 'urlopen', side_effect=lambda *a, **k: BytesIO(self.data)):
            fetch.prepare(self.cache, self.root)
        (self.cache / (self.names[0] + '.npz')).write_bytes(b'not trusted')
        with patch.object(fetch, 'urlopen', side_effect=AssertionError('unexpected network')):
            fetch.prepare(self.cache, self.root)
        with fetch.np.load(self.cache / (self.names[0] + '.npz')) as archive:
            self.assertEqual(archive['label'].tolist(), [1])
        self.assertEqual((self.root / 'catalog_inputs.csv').read_bytes(), self.receipt)

    def test_mismatch_does_not_touch_receipt_or_derived_products(self):
        for name in self.names:
            (self.cache / (name + '.npz')).write_bytes(b'old derived bytes')
        (self.cache / (self.names[0] + '.cat.ecsv')).write_bytes(self.data)
        (self.cache / (self.names[1] + '.cat.ecsv')).write_bytes(b'x' * len(self.data))
        with self.assertRaisesRegex(ValueError, 'byte mismatch'):
            fetch.prepare(self.cache, self.root)
        self.assertEqual((self.root / 'catalog_inputs.csv').read_bytes(), self.receipt)
        for name in self.names:
            self.assertEqual((self.cache / (name + '.npz')).read_bytes(), b'old derived bytes')

    def test_changed_frozen_url_is_rejected(self):
        receipt = self.root / 'catalog_inputs.csv'
        receipt.write_bytes(self.receipt.replace(b'https://stpubdata.s3.amazonaws.com/', b'https://example.invalid/'))
        with self.assertRaisesRegex(ValueError, 'invalid frozen catalog receipt'):
            fetch.prepare(self.cache, self.root)
        self.assertFalse(list(self.cache.iterdir()))

    def test_download_mismatch_and_wrong_membership_fail_before_outputs(self):
        with patch.object(fetch, 'urlopen', side_effect=lambda *a, **k: BytesIO(b'wrong')):
            with self.assertRaisesRegex(ValueError, 'byte mismatch'):
                fetch.prepare(self.cache, self.root)
        self.assertFalse(list(self.cache.iterdir()))
        self.assertEqual((self.root / 'catalog_inputs.csv').read_bytes(), self.receipt)
        with (self.root / 'pairs.csv').open('a') as f:
            f.write(self.names[0] + ',jw01208-o002_t001_nircam_clear-f200w\n')
        with self.assertRaisesRegex(ValueError, 'exactly match'):
            fetch.prepare(self.cache, self.root)
        self.assertEqual((self.root / 'catalog_inputs.csv').read_bytes(), self.receipt)


if __name__ == '__main__':
    unittest.main()
