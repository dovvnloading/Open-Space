import numpy as np,json,hashlib
from pathlib import Path
from astropy.io import fits
D=Path('D:/Open-Space/dossiers/gbm-qpo-20261117');angles=[(45.89,20.58),(45.11,45.31),(58.44,90.21),(314.87,45.24),(303.15,90.27),(3.35,89.79),(224.93,20.43),(224.62,46.18),(236.61,89.97),(135.19,45.55),(123.73,90.42),(183.74,90.32)]
def v(ra,dec):
 r,d=np.deg2rad([ra,dec]);return np.array([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])
out=[]
for p in sorted(D.glob('glg_trigdat*.fit')):
 with fits.open(p,memmap=True) as h:
  a=h[0].header;x=v(a['RA_SCX'],a['DEC_SCX']);z=v(a['RA_SCZ'],a['DEC_SCZ']);y=np.cross(z,x);source=v(a['RA_OBJ'],a['DEC_OBJ']);source_sc=np.array([source@x,source@y,source@z]);dets=[]
  for i,(az,ze) in enumerate(angles):
   az,ze=np.deg2rad([az,ze]);n=np.array([np.sin(ze)*np.cos(az),np.sin(ze)*np.sin(az),np.cos(ze)]);dets.append({'detector':i,'angle_deg':float(np.degrees(np.arccos(np.clip(n@source_sc,-1,1))))})
  ranked=sorted(dets,key=lambda q:q['angle_deg']);out.append({'trigger':p.name.split('_')[3],'source_radec':[a['RA_OBJ'],a['DEC_OBJ']],'axis_X':[a['RA_SCX'],a['DEC_SCX']],'axis_Z':[a['RA_SCZ'],a['DEC_SCZ']],'ranked':ranked,'note':'trigger-level pointing only, no time-dependent response or localization-error marginalization'});print(p.name,ranked[:3])
(D/'geometry-proposal.json').write_text(json.dumps(out,indent=2)+'\n')
