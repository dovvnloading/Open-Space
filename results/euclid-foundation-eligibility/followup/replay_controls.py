"""Recompute v5 grid-control scatter from deposited paired stamps."""
import argparse,numpy as np
p=argparse.ArgumentParser();p.add_argument('--stamps',required=True);a=p.parse_args();z=np.load(a.stamps);yy,xx=np.mgrid[-8:9,-8:9];disk=xx**2+yy**2<=16;ring=(xx**2+yy**2>=36)&(xx**2+yy**2<=64);allow=(1<<18)|(1<<24);v=[]
for n in range(len(z['blank_science1'])):
 m=[]
 for i in (1,2):
  im=z[f'blank_science{i}'][n];err=z[f'blank_rms{i}'][n];flag=z[f'blank_flag{i}'][n];good=np.isfinite(im)&np.isfinite(err)&(err>0)&((flag.astype(np.int64)&~allow)==0);assert good[disk].all() and good[ring].sum()>=50;m.append((float(np.sum(im[disk]-np.median(im[ring&good]))),float(np.sqrt(np.sum(err[disk]**2)))))
 v.append((m[1][0]-1.0061417496494092*m[0][0],float(np.hypot(m[1][1],1.0061417496494092*m[0][1]))))
x=np.array(v);print('count',len(x),'difference_sigma_mad',float(1.4826*np.median(np.abs(x[:,0]-np.median(x[:,0])))),'formal_median',float(np.median(x[:,1])))
