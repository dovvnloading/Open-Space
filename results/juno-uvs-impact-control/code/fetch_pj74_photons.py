import pathlib,json,time,subprocess,hashlib
p=pathlib.Path(__file__).resolve().parents[1]/'newer-pj74';u='https://pds-atmospheres.nmsu.edu/PDS/data/jnouvs_3001/DATA/ORBIT-74/UVS_S01_805506529_2025192_P74OBS_V01.FIT'
for _ in range(80):
 x=json.load(open(p/'hdu-layout.json'))
 if len(x)>2:break
 time.sleep(2)
h=x[2];assert h['width']==84 and h['rows']>=500000;a=h['data_offset'];b=a+500000*h['width']-1
subprocess.run(['curl','--fail','-Ls','--max-time','150','--retry','1','-r',f'{a}-{b}','-D',str(p/'pilot.http'),u,'-o',str(p/'pilot.bin')],check=True)
s=(p/'pilot.bin').read_bytes();assert len(s)==500000*h['width'];assert f'bytes {a}-{b}/' in (p/'pilot.http').read_text();json.dump({'source':u,'start':a,'end':b,'rows':500000,'sha256':hashlib.sha256(s).hexdigest()},open(p/'pilot-receipt.json','w'),indent=2);print('downloaded',len(s),flush=True)
