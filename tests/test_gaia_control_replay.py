"""Receipt-first control input replay regression tests."""
import hashlib
import importlib.util
from io import BytesIO
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / 'results/gaia-nss-dark-companions/orbit/fetch_controls.py'
spec = importlib.util.spec_from_file_location('gaia_control_fetch', MODULE)
fetch = importlib.util.module_from_spec(spec); spec.loader.exec_module(fetch)


class ControlReceipts(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name); self.cache=self.root/'cache';self.cache.mkdir();self.out=self.root/'derived';self.out.mkdir()
        self.urls={'one.dat':'https://example.test/one', 'two.dat':'https://example.test/two'}
        self.raw={'one.dat':b'first frozen input', 'two.dat':b'second frozen input'}
        self.receipts={n:dict(url=u,bytes=len(self.raw[n]),sha256=hashlib.sha256(self.raw[n]).hexdigest()) for n,u in self.urls.items()}
        self.receipt=self.root/'sources.json';self.receipt.write_text(json.dumps(self.receipts));self.before=self.receipt.read_bytes()
        for name in ['nss_controls.csv','published_rv.csv']:(self.out/name).write_bytes(b'original frozen table')
        self.patcher=patch.object(fetch,'EXPECTED_URLS',self.urls);self.patcher.start();self.addCleanup(self.patcher.stop)

    def parsed(self,*args):
        return ['source_id'],[{'source_id':'1'}],[['BH1','2450000','1','0.1','X']]

    def test_matching_cache_verifies_then_writes_only_separate_outputs(self):
        for n,b in self.raw.items():(self.cache/n).write_bytes(b)
        with patch.object(fetch,'parse_inputs',side_effect=self.parsed),patch.object(fetch,'urlopen',side_effect=AssertionError('network')):
            fetch.fetch(self.cache,self.out,self.receipt)
        self.assertIn(b'source_id',(self.out/'nss_controls.csv').read_bytes())
        self.assertIn(b'BH1',(self.out/'published_rv.csv').read_bytes())
        self.assertEqual(self.receipt.read_bytes(),self.before)

    def test_bad_cached_bytes_never_parse_or_rewrite_tables(self):
        (self.cache/'one.dat').write_bytes(self.raw['one.dat']);(self.cache/'two.dat').write_bytes(b'x'*len(self.raw['two.dat']))
        with patch.object(fetch,'parse_inputs',side_effect=AssertionError('parsed before verification')):
            with self.assertRaisesRegex(ValueError,'byte mismatch'):fetch.fetch(self.cache,self.out,self.receipt)
        self.assertEqual(self.receipt.read_bytes(),self.before)
        for n in ['nss_controls.csv','published_rv.csv']:self.assertEqual((self.out/n).read_bytes(),b'original frozen table')

    def test_bad_remote_bytes_do_not_write_cache_or_tables(self):
        with patch.object(fetch,'urlopen',side_effect=lambda url,**k:BytesIO(b'bad')):
            with self.assertRaisesRegex(ValueError,'byte mismatch'):fetch.fetch(self.cache,self.out,self.receipt)
        self.assertEqual(list(self.cache.iterdir()),[]);self.assertEqual(self.receipt.read_bytes(),self.before)
        for n in ['nss_controls.csv','published_rv.csv']:self.assertEqual((self.out/n).read_bytes(),b'original frozen table')

    def test_wrong_url_and_protected_source_output_are_rejected(self):
        self.receipts['one.dat']['url']='https://example.invalid/wrong'
        self.receipt.write_text(json.dumps(self.receipts))
        with self.assertRaisesRegex(ValueError,'invalid frozen receipt'):
            fetch.fetch(self.cache,self.out,self.receipt)
        for n in ['nss_controls.csv','published_rv.csv']:self.assertEqual((self.out/n).read_bytes(),b'original frozen table')
        with self.assertRaisesRegex(ValueError,'output must be separate'):
            fetch.fetch(self.cache,fetch.HERE,self.receipt)

    def test_parse_failure_and_bad_membership_leave_outputs_unchanged(self):
        for n,b in self.raw.items():(self.cache/n).write_bytes(b)
        with patch.object(fetch,'parse_inputs',side_effect=ValueError('bad table')):
            with self.assertRaisesRegex(ValueError,'bad table'):fetch.fetch(self.cache,self.out,self.receipt)
        for n in ['nss_controls.csv','published_rv.csv']:self.assertEqual((self.out/n).read_bytes(),b'original frozen table')
        self.receipts.pop('two.dat');self.receipt.write_text(json.dumps(self.receipts));changed=self.receipt.read_bytes()
        with self.assertRaisesRegex(ValueError,'membership mismatch'):fetch.fetch(self.cache,self.out,self.receipt)
        self.assertEqual(self.receipt.read_bytes(),changed)


if __name__=='__main__':unittest.main()
