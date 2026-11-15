"""Profile a smooth 2018 transition with a shared sigmoid and camera-specific amplitudes.
Every spatial component can vary in rate. The sigmoid width is the 10 to 90
percent rise time. The profile tests the sensitivity to contact morphology.
"""
from analyze import *
from scipy.optimize import minimize
from scipy.special import expit
from pathlib import Path
def term(v):
 n,p,g=read_visit(v); low,hi=58251.92,58252.43
 g=np.column_stack([np.maximum(g[:,0],low),np.minimum(g[:,1],hi)]);g=g[g[:,1]>g[:,0]]
 m=p.energy_ev.between(900,7000)&p.mjd_tdb_geo.between(low,hi)&(abs(p.x_arcsec-NEIGHBOR[0]/2)<19)&(abs(p.y_arcsec-NEIGHBOR[1]/2)<19)
 p=p[m];grid=(np.arange(90)+.5)*38/90-19;xx,yy=np.meshgrid(grid+NEIGHBOR[0]/2,grid+NEIGHBOR[1]/2)
 norm=components(xx.ravel(),yy.ravel(),4.5,1.5,19).mean(axis=1)*38**2
 X=components(p.x_arcsec.to_numpy(),p.y_arcsec.to_numpy(),4.5,1.5,19)/norm[:,None];return p.mjd_tdb_geo.to_numpy(),X,g
data=[term(f) for f in sorted((HERE/'inputs').glob('P0824450901*-photons.csv.gz'))]
def profile(par):
 mid=58252.+par[0]/100; s=(par[1]/100)/4.394449
 N=0
 for t,X,g in data:
  q=expit((t-mid)/s)
  h=(np.logaddexp(0,(g[:,1]-mid)/s)-np.logaddexp(0,(g[:,0]-mid)/s)).sum()*s*86.4
  ex=(g[:,1]-g[:,0]).sum()*86.4
  Z=np.vstack([X*(1-q),X*q]); dur=np.r_[np.ones(6)*(ex-h),np.ones(6)*h]
  def fun(b):
   y=np.maximum(b@Z,1e-12)
   return np.sum(b*dur)-np.log(y).sum(),dur-(Z/y).sum(axis=1)
  opt=minimize(fun,np.ones(12)*3,jac=True,method='L-BFGS-B',bounds=[(0,None)]*12,options={'ftol':1e-12});N+=opt.fun
 return N
if __name__=='__main__':
 op=minimize(profile,[13.45,3],method='Nelder-Mead',bounds=[(5,20),(.03,15)],options={'xatol':.005,'fatol':.0005})
 from scipy.optimize import minimize_scalar
 records=[]
 for width in [.03,.25,.5,.75,1.,1.5,2.,3.,4.,5.,7.,9.,12.,15.]:
  opt=minimize_scalar(lambda m:profile([m,width]),bounds=(8,18),method='bounded',options={'xatol':.008})
  records.append(dict(profile='width',mid_mjd=58252.+opt.x/100,width_days=width/100,delta_c=2*(opt.fun-op.fun)))
 for mid in np.arange(12,16.4,.2):
  opt=minimize_scalar(lambda w:profile([mid,w]),bounds=(2,16),method='bounded',options={'xatol':.06})
  records.append(dict(profile='midpoint',mid_mjd=58252.+mid/100,width_days=opt.x/100,delta_c=2*(opt.fun-op.fun)))
 records.append(dict(profile='optimum',mid_mjd=58252.+op.x[0]/100,width_days=op.x[1]/100,delta_c=0.))
 pd.DataFrame(records).to_csv(HERE/'xmm-ramp-profile.csv',index=False,float_format='%.8f')
