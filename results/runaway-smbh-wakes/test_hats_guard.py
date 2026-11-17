import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'catalog'))
from hats_screen import inspect
class CacheTests(unittest.TestCase):
 def test_existing_partition_rejected_even_matching_etag(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'12.json';p.write_text('{"etag":"same"}')
   with self.assertRaises(FileExistsError):inspect(('prefix/Npix=12',4,'same'),Path(tmp))
if __name__=='__main__':unittest.main()
