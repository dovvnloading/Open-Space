import concurrent.futures,subprocess,pathlib,io,math,json,re
import numpy as np
from astropy.io import fits
from astropy.time import Time
p=pathlib.Path(__file__).resolve().parents[1]/'clean-screen'
files=['UVS_S02_805506529_2025192_P74OBS_V01.FIT','UVS_S03_805506529_2025192_P74OBS_V01.FIT']+[f'UVS_S{i:02}_805533571_2025192_P74SY0_V01.FIT' for i in range(1,4)]
def run(name):
 u='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/'+name;f=p/(name+'.tail');hdr=p/(name+'.http')
 subprocess.run(['curl','--fail','-Ls','--retry','1','--max-time','90','-r','-12000000','-D',str(hdr),u,'-o',str(f)],check=True)
 b=f.read_bytes();start,end,total=map(int,re.findall(r'content-range: bytes (\d+)-(\d+)/(\d+)',hdr.read_text(),re.I)[-1]);assert len(b)==end-start+1
 out={'file':name,'bytes':len(b),'total_bytes':total,'hdus':[]}
 for match in re.finditer(b'XTENSION=',b):
  j=match.start()
  if (start+j)%2880:continue
  try:
   e=next(k+80 for k in range(j,len(b)-80,80) if b[k:k+8]==b'END     ');hs=math.ceil((e-j)/2880)*2880;h=fits.Header.fromstring(b[j:j+hs].decode(),sep='');name2=h.get('EXTNAME','');out['hdus'].append(name2)
   if name2!='Housekeeping Data':continue
   size=h['NAXIS1']*h['NAXIS2'];raw=b[j:j+hs+math.ceil(size/2880)*2880];data=fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode()+raw))[1].data
   np.savez_compressed(p/(name+'.hk.npz'),**{k:data[k] for k in data.names})
   et=data['EPHEMERIS_TIME'];clean=(data['SCI_OVERFLOW']==0)&(data['CHKSUM_ERROR']==0)&(data['HVPS_LIMITED']==0)&(data['HVPS_STATUS']==1)&(data['SAFETY_ACTIVE']==0)
   out.update({'hk_rows':len(et),'clean_housekeeping_rows':int(clean.sum()),'clean_indices':np.flatnonzero(clean).tolist(),'overflow_count':int((data['SCI_OVERFLOW']!=0).sum()),'stim_values':np.unique(data['DETECTOR_STIM']).tolist(),'utc_first':Time(2451545+et[0]/86400,format='jd',scale='tdb').utc.isot,'utc_last':Time(2451545+et[-1]/86400,format='jd',scale='tdb').utc.isot})
  except Exception as exc:out.setdefault('errors',[]).append(str(exc))
 json.dump(out,open(p/(name+'.summary.json'),'w'),indent=2);return out
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
 for result in ex.map(run,files):print(json.dumps(result),flush=True)
