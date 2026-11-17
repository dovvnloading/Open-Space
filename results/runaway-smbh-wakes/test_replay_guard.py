"""Reversible guard checks using tiny synthetic receipt and cache files."""
import csv,hashlib,json,tempfile,unittest
from pathlib import Path
from replay_guard import fresh_dirs,fresh_output,frozen_tiles,seal,sha,validate_cutouts,verify_seal,ROOT

class GuardTests(unittest.TestCase):
 def test_input_mutation(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);rows=[]
   for kind in ('VIS','rms'):
    p=root/('1_'+kind+'.fits');p.write_bytes(kind.encode());rows.append({'tile':'1','product':kind,'bytes':p.stat().st_size,'sha256':sha(p)})
   manifest=root/'receipts.csv'
   with manifest.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
   self.assertTrue(validate_cutouts(root,['1'],manifest));(root/'1_VIS.fits').write_bytes(b'bad')
   with self.assertRaises(ValueError):validate_cutouts(root,['1'],manifest)
 def test_frozen_selection_and_change(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'tiles.csv';p.write_bytes((ROOT/'euclid/q1_tiles.csv').read_bytes());self.assertEqual(len(frozen_tiles(p)),344);p.write_bytes(p.read_bytes()+b'x')
   with self.assertRaises(ValueError):frozen_tiles(p)
 def test_collision_normalization_and_checkpoint(self):
  with tempfile.TemporaryDirectory() as tmp:
   base=Path(tmp);source=base/'input';source.mkdir()
   with self.assertRaises(ValueError):fresh_output(source,source/'child')
   with self.assertRaises(ValueError):fresh_dirs(base/'same',base/'a'/'..'/'same')
   dest=fresh_output(source,base/'new');(dest/'one.json').write_bytes(b'one');seal(dest,'screen');self.assertEqual(len(verify_seal(dest,'screen')),1)
   with self.assertRaises(FileExistsError):fresh_dirs(dest)
   manifest=dest/'checkpoint_manifest.json';saved=manifest.read_bytes();record=json.loads(saved);record['code_sha256']['replay_guard.py']='0'*64;manifest.write_text(json.dumps(record))
   with self.assertRaises(ValueError):verify_seal(dest,'screen')
   manifest.write_bytes(saved);(dest/'one.json').write_bytes(b'two')
   with self.assertRaises(ValueError):verify_seal(dest,'screen')
   with self.assertRaises(ValueError):absolute_relative()

def absolute_relative():
 from replay_guard import absolute_external
 return absolute_external(Path('relative'))
if __name__=='__main__':unittest.main()
