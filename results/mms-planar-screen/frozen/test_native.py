import tempfile,unittest,sys,json
from pathlib import Path
from unittest.mock import patch
import cdflib,numpy as np
import development
from plane import proposal

def write(path,variables):
 w=cdflib.cdfwrite.CDF(path,delete=True)
 for name,data,units in variables:
  data=np.asarray(data);dims=list(data.shape[1:]);kind=33 if name=='Epoch' else (4 if 'flag' in name else 45)
  attrs={'UNITS':units,'FILLVAL':-1e31,'VALIDMIN':-20000.,'VALIDMAX':20000.} if kind==45 else {'UNITS':units}
  if 'err_brst' in name:attrs.update(VALIDMIN=0.,VALIDMAX=5000.)
  if 'bulkv_gse' in name:attrs.update(VALIDMIN=-2500.,VALIDMAX=2500.)
  if 'fgm_b_' in name:attrs.update(VALIDMIN=[-20000.,-20000.,-20000.,0.],VALIDMAX=[20000.]*4)
  w.write_var({'Variable':name,'Data_Type':kind,'Num_Elements':1,'Rec_Vary':True,'Dim_Sizes':dims},var_attrs=attrs,var_data=data)
 w.close()
class NativeTests(unittest.TestCase):
 def test_native_interval_gaps_flags_fill(self):
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as tmp:
   path=Path(tmp);start=cdflib.cdfepoch.compute_tt2000([2026,9,30,13,28,0,0,0,0]);tt=start+np.arange(2400,dtype=np.int64)*10000000;dt=start+np.arange(160,dtype=np.int64)*150000000
   t=np.arange(2400)/100.;a=np.clip(t-9.5,0,1)*np.pi/2;b=np.column_stack([8*np.cos(a),8*np.sin(a),.02*np.sin(t*5),np.ones(len(t))*8]);b[500,0]=-1e31
   ff=np.zeros(len(t),np.int32);ff[600]=1;dv=np.tile([0.,0.,100.],(len(dt),1));de=np.ones_like(dv);flags=np.zeros(len(dt),np.int32);flags[50]=1;flags[60]=1
   dt=np.delete(dt,70);dv=np.delete(dv,70,axis=0);de=np.delete(de,70,axis=0);flags=np.delete(flags,70)
   write(path/'f.cdf',[('Epoch',tt,'ns'),('mms1_fgm_b_gse_brst_l2',b,'nT'),('mms1_fgm_flag_brst_l2',ff,' '),('mms1_fgm_r_gse_brst_l2',np.zeros((len(tt),3)),'km')])
   write(path/'d.cdf',[('Epoch',dt,'ns'),('mms1_dis_bulkv_gse_brst',dv,'km/s'),('mms1_dis_bulkv_err_brst',de,'km/s'),('mms1_dis_errorflags_brst',flags,' ')])
   f=cdflib.CDF(path/'f.cdf');d=cdflib.CDF(path/'d.cdf');origin,ts,bb,v,e,good=development.read_sc1(f,d,development.EXPECTED)
   self.assertGreaterEqual(origin,start+3413441756);self.assertLessEqual(origin+round(ts[-1]*1e9),start+24157399823);self.assertTrue((~good).any());self.assertTrue(good.any());self.assertEqual(bb.shape[1],3);self.assertTrue(np.all(np.diff(ts)>0))
   self.assertIsNone(proposal(ts,bb,v,good));del f,d
   import gc;gc.collect()
 def test_auth_before_science(self):
  argv=['development.py','--manifest','unused','--manifest-sha256','m','--authorization','unused','--authorization-sha256','bad','--output','D:/unused']
  with patch.object(sys,'argv',argv),patch.object(development,'verify',return_value={}),patch.object(development,'authorize',side_effect=ValueError('sealed')),patch.object(development.cdflib,'CDF',side_effect=AssertionError('science opened')):
   with self.assertRaisesRegex(ValueError,'sealed'):development.main()
 def test_source_only_positions_guard(self):
  class Source:
   def __init__(self,craft):self.craft=craft
   def varget(self,key):
    if self.craft!=1 and key not in ('Epoch_state',f'mms{self.craft}_fgm_r_gse_brst_l2'):raise AssertionError('heldout physical value accessed')
    if key=='Epoch_state':return np.array([0,10000000000],dtype=np.int64)
    return np.zeros((2,3))
  class Dirs:
   def __truediv__(self,n):return Path(n)
  ids={f'mms{i}_fgm_test.cdf':'x' for i in range(1,5)};ids['mms1_fpi_test.cdf']='x'
  t=np.arange(0,10,.01);v=np.ones((len(t),3));x={'index':5,'time':5.,'normal':[1,0,0],'normal_speed_km_s':50.}
  with tempfile.TemporaryDirectory(dir='D:/Open-Space/campaigns') as folder:
   args=['development.py','--manifest','unused','--manifest-sha256','m','--authorization','unused','--authorization-sha256','a','--output',str(Path(folder)/'out.json')]
   def getcdf(path):
    name=str(path);return Source(int(name[3])) if name.startswith('mms') else Source(1)
   with patch.object(sys,'argv',args),patch.object(development,'verify',return_value={'inputs':ids,'git_head':'synthetic'}),patch.object(development,'authorize'),patch.object(development,'RAW',Dirs()),patch.object(development.cdflib,'CDF',side_effect=getcdf),patch.object(development,'read_sc1',side_effect=lambda f,d,c,a:(a.update(cropped_first_fgm_native_index=0) or (0,t,v,v,v,np.ones(len(t),bool)))),patch.object(development,'proposal',side_effect=lambda t,b,v,g,a:(a.update(quality_cadence_windows=1,rotation_windows=1) or x)),patch.object(development,'sensitivity',return_value={'status':'inconclusive_normal_speed'}):development.main()
   self.assertEqual(json.loads((Path(folder)/'out.json').read_text())['status'],'inconclusive_normal_speed')
if __name__=='__main__':unittest.main()
