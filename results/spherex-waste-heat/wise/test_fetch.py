"""Small immutable-cache regressions without network access."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('wise_fetch',Path(__file__).with_name('fetch.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
T={'label':'A','ra':'10','dec':'20','ref_epoch':'2016','pmra':'1','pmdec':'2'}
class Response:
    status=200
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def read(self):return (m.TABLES['allwise'][1]+'\n').encode()
class CacheTest(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.p=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def create(self):
        with patch.object(m,'urlopen',return_value=Response()):m.fetch('allwise',T,self.p)
    def test_zero_rows_and_matching_pair(self):
        self.create();self.assertEqual(json.loads((self.p/'A_allwise.json').read_text())['rows'],0)
        with patch.object(m,'urlopen',side_effect=AssertionError('network used')):m.fetch('allwise',T,self.p)
    def test_changed_query(self):
        self.create();before={p.name:p.read_bytes() for p in self.p.iterdir()}
        with self.assertRaises(ValueError):m.fetch('allwise',dict(T,ra='11'),self.p)
        self.assertEqual(before,{p.name:p.read_bytes() for p in self.p.iterdir()})
    def test_hash_or_size_mismatch(self):
        self.create();p=self.p/'A_allwise.json';r=json.loads(p.read_text());r['bytes']+=1;p.write_text(json.dumps(r));before={x.name:x.read_bytes() for x in self.p.iterdir()}
        with self.assertRaises(ValueError):m.fetch('allwise',T,self.p)
        self.assertEqual(before,{x.name:x.read_bytes() for x in self.p.iterdir()})
    def test_wrong_digest(self):
        self.create();p=self.p/'A_allwise.json';r=json.loads(p.read_text());r['sha256']='0'*64;p.write_text(json.dumps(r));before={x.name:x.read_bytes() for x in self.p.iterdir()}
        with self.assertRaises(ValueError):m.fetch('allwise',T,self.p)
        self.assertEqual(before,{x.name:x.read_bytes() for x in self.p.iterdir()})
    def test_partial_pair(self):
        self.create();(self.p/'A_allwise.json').unlink();before=(self.p/'A_allwise.csv.gz').read_bytes()
        with self.assertRaises(ValueError):m.fetch('allwise',T,self.p)
        self.assertEqual(before,(self.p/'A_allwise.csv.gz').read_bytes())
    def test_error_preserved(self):
        with patch.object(m,'urlopen',side_effect=OSError('offline')):
            with self.assertRaises(OSError):m.fetch('allwise',T,self.p)
        p=self.p/'A_allwise.error.json';before=p.read_bytes()
        with self.assertRaises(ValueError):m.fetch('allwise',T,self.p)
        self.assertEqual(before,p.read_bytes())
    def test_fetch_to_measure_workflow(self):
        import sys
        sys.modules['fetch']=m
        spec2=importlib.util.spec_from_file_location('wise_measure',Path(__file__).with_name('measure.py'))
        measure=importlib.util.module_from_spec(spec2);spec2.loader.exec_module(measure)
        archive=self.p/'archive';out=self.p/'derived'
        class BandResponse(Response):
            def __init__(self,kind):self.kind=kind
            def read(self):return (m.TABLES[self.kind][1]+'\n').encode()
        for kind in m.TABLES:
            with patch.object(m,'urlopen',return_value=BandResponse(kind)):
                m.fetch(kind,T,archive)
        before={x.name:x.read_bytes() for x in archive.iterdir()}
        with patch.object(measure,'targets',return_value=[dict(T,source_id='123')]):
            measure.main(archive,out)
        self.assertEqual(before,{x.name:x.read_bytes() for x in archive.iterdir()})
        self.assertEqual({x.name for x in out.iterdir()},{'catalog_photometry.csv','visits.csv','variability.csv','queries.json'})
        with self.assertRaises(FileExistsError):measure.main(archive,out)
if __name__=='__main__':unittest.main()
