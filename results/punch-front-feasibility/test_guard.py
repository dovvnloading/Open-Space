import copy, json, unittest
from pathlib import Path
from guard import verify
ROOT=Path(__file__).resolve().parent
class GuardTests(unittest.TestCase):
 def setUp(self):
  self.prov=json.loads((ROOT/'provenance.json').read_text()); self.parents=json.loads((ROOT/'parent_headers.json').read_text())
 def test_metadata(self): self.assertEqual(verify(self.prov,self.parents)['parents'],9)
 def test_science_denied(self):
  with self.assertRaises(ValueError): verify(self.prov,self.parents,'development')
 def test_parent_reuse_denied(self):
  self.prov[1]['parents'][0]=self.prov[0]['parents'][0]
  with self.assertRaises(ValueError): verify(self.prov,self.parents)
 def test_time_overlap_denied(self):
  self.parents[3]['headers']['DATE-BEG']='2026-11-01T00:00:20'
  with self.assertRaises(ValueError): verify(self.prov,self.parents)
 def test_wrong_instrument_denied(self):
  self.parents[1]['headers']['OBSCODE']='4'
  with self.assertRaises(ValueError): verify(self.prov,self.parents)
if __name__=='__main__': unittest.main()
