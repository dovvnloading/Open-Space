import pathlib,subprocess,sys,tempfile,shutil,unittest,os,importlib.util,io
ROOT=pathlib.Path(__file__).resolve().parents[1]/'results/fast-extragalactic-alerts-2026'
SCRIPT=ROOT/'rubin/screen.py'
class FrozenRubin(unittest.TestCase):
 def test_frozen_input_output_and_units(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp); data=p/'data';shutil.copytree(ROOT/'rubin/data',data)
   def run(out,*extra):return subprocess.run([sys.executable,str(SCRIPT),str(data),str(out),*extra],capture_output=True)
   out=p/'ok';self.assertEqual(run(out).returncode,0); before=(out/'funnel.json').read_bytes();self.assertNotEqual(run(out).returncode,0);self.assertEqual(before,(out/'funnel.json').read_bytes())
   self.assertNotEqual(run(p/'badunit','--flux-unit','Jy').returncode,0);self.assertFalse((p/'badunit').exists())
   f=data/'day_61349.csv.gz'; original=f.read_bytes();f.write_bytes(original+b'bad');self.assertNotEqual(run(p/'corrupt').returncode,0);self.assertFalse((p/'corrupt').exists());f.write_bytes(original)
   (data/'day_61350.csv.gz').unlink();self.assertNotEqual(run(p/'missing').returncode,0);self.assertFalse((p/'missing').exists())
 def test_ztf_rejects_modified_receipt_inputs_before_output(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=pathlib.Path(tmp); z=p/'ztf';shutil.copytree(ROOT/'ztf',z)
   script=z/'scripts/ztf_screen.py'; out=p/'result'
   def run():return subprocess.run([sys.executable,str(script),str(out)],capture_output=True)
   self.assertEqual(run().returncode,0); before=(out/'ztf_counts.json').read_bytes();self.assertNotEqual(run().returncode,0);self.assertEqual(before,(out/'ztf_counts.json').read_bytes())
   shutil.rmtree(out); f=z/'data/ztf/61349_detections.csv.gz';f.write_bytes(f.read_bytes()+b'bad');self.assertNotEqual(run().returncode,0);self.assertFalse(out.exists())
 def test_normalized_live_campaign_and_existing_response(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=pathlib.Path(tmp)/'campaign';(root/'data/ztf').mkdir(parents=True);(root/'receipts/ztf').mkdir(parents=True)
   old=os.environ.get('FAST_ALERTS_NEW_CAMPAIGN');os.environ['FAST_ALERTS_NEW_CAMPAIGN']=str(root/'..'/'campaign')
   try:
    spec=importlib.util.spec_from_file_location('tap_test',ROOT/'ztf/scripts/ztf_tap.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    class Reply(io.BytesIO):
     status=200;headers={'Content-Type':'text/csv'}
    calls=[]
    def fake(*args,**kwargs):calls.append(1);return Reply(b'nrows\n0\n')
    mod.urllib.request.urlopen=fake
    mod.get('check','SELECT 1',mod.RAW);self.assertEqual(len(calls),1)
    with self.assertRaises(FileExistsError):mod.get('check','SELECT 1',mod.RAW)
    self.assertEqual(len(calls),1)
   finally:
    if old is None:os.environ.pop('FAST_ALERTS_NEW_CAMPAIGN',None)
    else:os.environ['FAST_ALERTS_NEW_CAMPAIGN']=old
if __name__=='__main__':unittest.main()
