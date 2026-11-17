import pymupdf,pathlib,re,csv,numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
p=pathlib.Path('D:/Open-Space/campaigns/vast-tde-expanded-20261117');d=pymupdf.open(p/'vast-paper.pdf');a=[]
for page in [13,14,15]:a+=re.findall(r'ASKAP\s+J(\d{6}(?:\.\d)?)([+−-])(\d{6})',d[page].get_text())
print('names',len(a),'unique',len(set(a)));(p/'vast-table2-names.txt').write_text('\n'.join('ASKAP J'+h+sign+dec for h,sign,dec in a)+'\n',encoding='utf-8')
coords=SkyCoord([f'{h[:2]}:{h[2:4]}:{h[4:]}' for h,s,d in a],[f'{"-" if s in "−-" else "+"}{d[:2]}:{d[2:4]}:{d[4:]}' for h,s,d in a],unit=(u.hourangle,u.deg))
for name,ra,dec in [('AT2020nov','16:58:12.970','+02:07:03.04'),('AT2018hyz','10:06:50.871','+01:41:34.08'),('ASASSN-15oi','20:39:09.096','-30:45:20.71')]:
 t=SkyCoord(ra,dec,unit=(u.hourangle,u.deg));idx,sep,_=t.match_to_catalog_sky(coords);print(name,'nearest table2',a[int(np.asarray(idx).ravel()[0])],float(np.asarray(sep.arcsec).ravel()[0]),'arcsec')
