from pathlib import Path
from astropy.io import fits
from scipy.ndimage import maximum_filter
import numpy as np
p=Path(r'D:\Open-Space\campaigns\euclid-foundation-20261116-70b26a3-terra61\cutouts'); base='102158276_013663548_VIS-'
a=fits.getdata(p/(base+'science.fits')); e=fits.getdata(p/(base+'rms.fits')); f=fits.getdata(p/(base+'flag.fits')); good=np.isfinite(a)&np.isfinite(e)&(e>0)&(f==0); peaks=np.argwhere(good & (a>8*e) & (a==maximum_filter(a,size=3))); selected=[]
for y,x in sorted(peaks,key=lambda v:float(a[tuple(v)]),reverse=True):
 if min(x,y,299-x,299-y)<10:continue
 yy,xx=np.ogrid[-4:5,-4:5];disk=xx**2+yy**2<=16
 if not good[y-4:y+5,x-4:x+5][disk].all():continue
 if any((x-x1)**2+(y-y1)**2<100 for x1,y1 in selected):continue
 selected.append((int(x),int(y)))
print('Eligible maxima:',len(selected));print(selected)
