import subprocess,pathlib,math,json
from astropy.io import fits
p=pathlib.Path(__file__).resolve().parents[1]/'newer-pj74';p.mkdir(exist_ok=True)
u='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/UVS_S01_805506529_2025192_P74OBS_V01.FIT'
o=0;infos=[]
for i in range(9):
 f=p/f'header-{i}.bin';http=p/f'header-{i}.http'
 subprocess.run(['curl','--fail','-Ls','--retry','1','--max-time','40','-r',f'{o}-{o+287999}','-D',str(http),u,'-o',str(f)],check=True)
 b=f.read_bytes();assert len(b)>2880;e=next(j+80 for j in range(0,len(b),80) if b[j:j+8]==b'END     ');hs=math.ceil(e/2880)*2880;h=fits.Header.fromstring(b[:hs].decode(),sep='');(p/f'header-{i}.txt').write_text(h.tostring());n=abs(h['BITPIX'])//8
 for j in range(h['NAXIS']):n*=h[f'NAXIS{j+1}']
 n+=h.get('PCOUNT',0);info={'hdu':i,'offset':o,'header_size':hs,'data_offset':o+hs,'data_bytes':n,'name':h.get('EXTNAME','primary'),'rows':h.get('NAXIS2'),'width':h.get('NAXIS1')};infos.append(info);print(info,flush=True);json.dump(infos,open(p/'hdu-layout.json','w'),indent=2);o+=hs+math.ceil(n/2880)*2880
