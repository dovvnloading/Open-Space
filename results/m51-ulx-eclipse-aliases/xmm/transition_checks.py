"""Test assignment of the 2018 change to ULX-1 with free companion/background."""
from analyze import *
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2
def check(visit,mid,lo=900,hi=7000,shift=0,rc=4.5,alpha=1.5):
    name,p,gti=visit
    start,stop=58251.92,58252.43
    gti=np.column_stack((np.maximum(gti[:,0],start),np.minimum(gti[:,1],stop)))
    gti=gti[gti[:,1]>gti[:,0]]
    disp=shift*NEIGHBOR/np.linalg.norm(NEIGHBOR)
    xx,yy=p.x_arcsec-disp[0],p.y_arcsec-disp[1]
    take=p.energy_ev.between(lo,hi)&p.mjd_tdb_geo.between(start,stop)&(abs(xx-NEIGHBOR[0]/2)<19)&(abs(yy-NEIGHBOR[1]/2)<19)
    x=xx[take].to_numpy();y=yy[take].to_numpy();t=p.mjd_tdb_geo[take].to_numpy()
    grid=(np.arange(90)+.5)*38/90-19;gx,gy=np.meshgrid(grid+NEIGHBOR[0]/2,grid+NEIGHBOR[1]/2)
    norm=components(gx.ravel(),gy.ravel(),rc,alpha,19).mean(axis=1)*38**2
    X=components(x,y,rc,alpha,19)/norm[:,None]
    pre=t<mid; dur=np.array([np.maximum(0,np.minimum(gti[:,1],mid)-gti[:,0]).sum()*86.4,0]);dur[1]=(gti[:,1]-gti[:,0]).sum()*86.4-dur[0]
    nn=[]; errors=[]; null=[]
    for a,w in enumerate([pre,~pre]):
        raw=fit(x[w],y[w],rc,alpha)
        nn.append(np.array(raw[:3])/dur[a])
    def solve(tied):
        npar=11 if tied else 12
        def fun(a):
            if tied:
                pars=[a[:6],np.r_[a[0],a[6:]]]
            else:
                pars=[a[:6],a[6:]]
            N=0;grad=[]
            for k,m in enumerate([pre,~pre]):
                model=np.maximum(pars[k]@X[:,m],1e-12)
                N+=dur[k]*pars[k].sum()-np.log(model).sum()
                grad.append(dur[k]-(X[:,m]/model).sum(axis=1))
            if tied:
                grad=np.r_[grad[0][0]+grad[1][0],grad[0][1:],grad[1][1:]]
            else:grad=np.r_[grad[0],grad[1]]
            return N,grad
        res=minimize(fun,np.full(npar,4.),jac=True,method='L-BFGS-B',bounds=[(0,None)]*npar,options={'ftol':1e-12})
        return res.fun
    ll0,ll1=solve(True),solve(False)
    a,b=nn
    return dict(detector=name[11:13],pre_ks=dur[0],post_ks=dur[1],ulx1_pre_rate=a[0],ulx1_pre_error=a[1],ulx1_post_rate=b[0],ulx1_post_error=b[1],ulx2_pre_rate=a[2],ulx2_post_rate=b[2],delta_c_vary_ulx1=2*(ll0-ll1))
if __name__=='__main__':
    mid=58252.1345; rows=[]
    for ff in sorted((HERE/'inputs').glob('P0824450901*-photons.csv.gz')):
        v=read_visit(ff)
        for label,kw in [('base',{}),('soft',dict(lo=500,hi=1500)),('hard',dict(lo=1500,hi=7000)),('shift_plus',dict(shift=1.5)),('shift_minus',dict(shift=-1.5)),('broad',dict(rc=6.,alpha=1.4)),('narrow',dict(rc=3.5,alpha=1.65))]:
            rows.append(dict(variation=label,**check(v,mid,**kw)))
    f=pd.DataFrame(rows)
    f.to_csv(HERE/'xmm-transition-checks.csv',index=False,float_format='%.7g')
    print(f[f.variation=='base'].round(2).to_string(index=False))
