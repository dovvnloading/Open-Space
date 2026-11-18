import hashlib,json,os,pathlib,subprocess,sys,tempfile,unittest
REPO=pathlib.Path(__file__).resolve().parents[1]
RESULT=REPO/'results/stellar-stokes-v'
SCRIPT=RESULT/'validation/dispositions.py'
@unittest.skipUnless(os.name == 'nt' and pathlib.Path('D:/').is_dir(), 'requires Windows D: replay destination')
class ReplayIsolation(unittest.TestCase):
 def test_external_snapshot_and_collision(self):
  with tempfile.TemporaryDirectory() as tmp:
   out=pathlib.Path(tmp)/'run'
   command=[sys.executable,str(SCRIPT),'--out',str(out)]
   originals={n:(RESULT/'validation'/n).read_bytes() for n in ('match_dispositions.csv','host_v_measurements.csv','disposition_receipt.json')}
   first=subprocess.run(command,capture_output=True,text=True);self.assertEqual(first.returncode,0,first.stderr)
   for n,b in originals.items():self.assertEqual((out/n).read_bytes(),b);self.assertEqual((RESULT/'validation'/n).read_bytes(),b)
   manifest=json.loads((out/'input_manifest.json').read_text());self.assertEqual(len(manifest['inputs']),15)
   self.assertEqual(manifest['code_sha256'],hashlib.sha256(SCRIPT.read_bytes()).hexdigest())
   for n,info in manifest['inputs'].items():self.assertEqual(hashlib.sha256((out/'inputs'/n).read_bytes()).hexdigest(),info['sha256'])
   before=(out/'input_manifest.json').read_bytes();self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0);self.assertEqual((out/'input_manifest.json').read_bytes(),before)
   self.assertNotEqual(subprocess.run([sys.executable,str(SCRIPT),'--out',str(RESULT/'unsafe')],capture_output=True).returncode,0)
   self.assertNotEqual(subprocess.run([sys.executable,str(SCRIPT),'--out','relative-output'],capture_output=True).returncode,0)
   self.assertNotEqual(subprocess.run([sys.executable,str(SCRIPT),'--out',r'C:\stokes-replay-forbidden'],capture_output=True).returncode,0)
if __name__=='__main__':unittest.main()
