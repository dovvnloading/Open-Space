import numpy as np,json,hashlib
from pathlib import Path
from astropy.io import fits
from scipy.spatial.transform import Rotation
D=Path('D:/Open-Space/dossiers/gbm-qpo-20261117');angs=[(45.89,20.58),(45.11,45.31),(58.44,90.21),(314.87,45.24),(303.15,90.27),(3.35,89.79),(224.93,20.43),(224.62,46.18),(236.61,89.97),(135.19,45.55),(123.73,90.42),(183.74,90.32)];spec={'261002867':(155.,-1.9,5.5,8.96016),'261007236':(11.61403,-17.69197,.00075,22.7844)};out=[]
for bid,(ra,dec,unc,duration) in spec.items():
 p=D/f'glg_trigdat_all_bn{bid}_v01.fit';r,d=np.deg2rad([ra,dec]);source=np.array([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])
 with fits.open(D/f'glg_bcat_all_bn{bid}_v00.fit',memmap=True) as bc: offset=float(bc[0].header['T90START'])
 with fits.open(p,memmap=True) as h:
  table=h['EVNTRATE'].data;trigger=h[0].header['TRIGTIME'];tt=np.asarray(table['TIME']);q=np.asarray(table['SCATTITD']);idx=np.where((tt>=trigger-1)&(tt<=trigger+offset+duration+1))[0];
  if len(idx)<2:raise ValueError('no attitude coverage')
  mats=Rotation.from_quat(q[idx]).as_matrix();v=np.einsum('nji,j->ni',mats,source);rank=[]
  for i,(az,ze) in enumerate(angs):
   aa,zz=np.deg2rad([az,ze]);normal=np.array([np.sin(zz)*np.cos(aa),np.sin(zz)*np.sin(aa),np.cos(zz)]);ar=np.rad2deg(np.arccos(np.clip(v@normal,-1,1)));rank.append({'detector':i,'angle_min':float(ar.min()),'angle_max':float(ar.max())})
  order=sorted(rank,key=lambda x:x['angle_max']);out.append({'trigger':bid,'source_radec':[ra,dec],'source_uncertainty_deg':unc,'attitude_rows':len(idx),'attitude_time_offsets_s':[float((tt[idx]-trigger).min()),float((tt[idx]-trigger).max())],'detector_ranked':order,'trigdat_sha256':hashlib.sha256(p.read_bytes()).hexdigest()});print(bid,'rank',order[:3],'attitude',len(idx))
(D/'refined-geometry.json').write_text(json.dumps(out,indent=2)+'\n')
