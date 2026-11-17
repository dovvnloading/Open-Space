import csv
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
FILE = ROOT / 'results/jwst-repeat-imaging-transients/goods-s-oct2026/acquire.py'
spec = importlib.util.spec_from_file_location('goods_acquire', FILE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class ReceiptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.receipt = self.root / 'receipts'; self.receipt.mkdir()
        self.cache = self.root / 'cache'
        self.key = 'jwst/public/jw08060/L3/t/o001/a_cat.ecsv'
        self.row = {'local_path':'8060/a_cat.ecsv', 'url':m.BASE+self.key, 'bytes':'3', 'sha256':hashlib.sha256(b'abc').hexdigest()}
        self.write()
    def tearDown(self): self.tmp.cleanup()
    def write(self, rows=None, keys=None):
        for name, fields, data in [('8060_manifest.csv',list(self.row),rows if rows is not None else [self.row]),('8060_index.csv',['key','bytes'],keys if keys is not None else [{'key':self.key,'bytes':'3'}])]:
            with (self.receipt/name).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(data)
    def run_acquire(self, fetch): return m.acquire('8060',self.cache,self.receipt,fetch)
    def test_matching_cache_offline(self):
        p=self.cache/'8060/a_cat.ecsv'; p.parent.mkdir(parents=True); p.write_bytes(b'abc')
        before=(self.receipt/'8060_manifest.csv').read_bytes()
        self.assertEqual(self.run_acquire(lambda _: self.fail('network called')),1)
        self.assertEqual((self.receipt/'8060_manifest.csv').read_bytes(),before)
    def test_bad_remote_or_truncated_cache_no_changes(self):
        before=(self.receipt/'8060_manifest.csv').read_bytes()
        with self.assertRaises(ValueError): self.run_acquire(lambda _: b'bad')
        self.assertFalse(self.cache.exists())
        p=self.cache/'8060/a_cat.ecsv'; p.parent.mkdir(parents=True); p.write_bytes(b'ab')
        with self.assertRaises(ValueError): self.run_acquire(lambda _: self.fail('network called'))
        self.assertEqual(p.read_bytes(),b'ab'); self.assertEqual((self.receipt/'8060_manifest.csv').read_bytes(),before)
    def test_duplicate_traversal_and_membership(self):
        for rows,keys in [([self.row,self.row],None),([{**self.row,'local_path':'8060/../a_cat.ecsv'}],None),([self.row],[{'key':self.key,'bytes':'4'}]),([] ,None)]:
            self.write(rows,keys)
            with self.assertRaises(ValueError): self.run_acquire(lambda _: b'abc')
            self.assertFalse(self.cache.exists())
    def test_all_validate_before_cache_write(self):
        key2=self.key.replace('a_cat','b_cat'); row2={**self.row,'local_path':'8060/b_cat.ecsv','url':m.BASE+key2}
        self.write([self.row,row2],[{'key':self.key,'bytes':'3'},{'key':key2,'bytes':'3'}])
        with self.assertRaises(ValueError): self.run_acquire(lambda url: b'abc' if url.endswith('a_cat.ecsv') else b'bad')
        self.assertFalse(self.cache.exists())

if __name__=='__main__': unittest.main()
