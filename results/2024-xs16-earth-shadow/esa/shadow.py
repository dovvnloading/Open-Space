from pathlib import Path
from datetime import datetime,timedelta
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.stats import qmc, norm

HERE=Path(__file__).resolve().parent
AU=149597870.7
GM=0.0002959122082855911
MASSES={1:6023600.,2:408523.71,399:332946.0487,301:27068703.,4:3098708.,5:1047.3486,6:3497.898,7:22902.98,8:19412.24}
EPOCH=60655.500491011  # MJD TT, mid-arc epoch of frozen Aegis solution
utc0=datetime(2024,12,11)
def tdb(s):
    return (s+69.184)/86400

def orbit():
    rows=(HERE/'2024XS16.eq0').read_text().splitlines()
    vals=np.array(next(l for l in rows if l.startswith(' EQU')).split()[1:],dtype=float)
    tri=np.fromstring(' '.join(' '.join(l.split()[1:]) for l in rows if l.startswith(' COV')),sep=' ')
    cov=np.zeros((6,6));cov[np.triu_indices(6)]=tri;cov+=np.triu(cov,1).T
    return vals,cov

def elements_to_state(rows):
    a,h,k,p,q,lam=rows.T
    e=np.hypot(h,k);varpi=np.arctan2(h,k);o=np.arctan2(p,q);w=varpi-o
    inc=2*np.arctan(np.hypot(p,q));m=(np.deg2rad(lam)-varpi+np.pi)%(2*np.pi)-np.pi
    E=m.copy()
    for _ in range(12):E-=(E-e*np.sin(E)-m)/(1-e*np.cos(E))
    x=a*(np.cos(E)-e);y=a*np.sqrt(1-e*e)*np.sin(E)
    n=np.sqrt(GM/a)/(1-e*np.cos(E));vx=-n*np.sin(E);vy=n*np.sqrt(1-e*e)*np.cos(E)
    c,s=np.cos,np.sin
    P=np.array([c(o)*c(w)-s(o)*s(w)*c(inc),s(o)*c(w)+c(o)*s(w)*c(inc),s(w)*s(inc)]).T
    Q=np.array([-c(o)*s(w)-s(o)*c(w)*c(inc),-s(o)*s(w)+c(o)*c(w)*c(inc),c(w)*s(inc)]).T
    return np.column_stack((x[:,None]*P+y[:,None]*Q,vx[:,None]*P+vy[:,None]*Q))

def planets():
    spl={}
    data=np.loadtxt(HERE/'planets.csv',delimiter=',',skiprows=1)
    for n in MASSES:
        rows=data[data[:,0]==n,1:]
        jd=rows[:,0]-2460655.5
        spl[n]=CubicHermiteSpline(jd,rows[:,1:4],rows[:,4:7],extrapolate=False)
    return spl

def integrate(elements,spl):
    y0=elements_to_state(elements)
    def force(t,y):
        r,v=np.split(y.reshape(-1,6),2,axis=1)
        acc=-GM*r/np.linalg.norm(r,axis=1)[:,None]**3
        for n,ratio in MASSES.items():
            p=spl[n](t);d=p-r
            acc+=GM/ratio*(d/np.linalg.norm(d,axis=1)[:,None]**3-p/np.linalg.norm(p)**3)
        return np.column_stack((v,acc)).ravel()
    sol=solve_ivp(force,(EPOCH-60655,.86),y0.ravel(),method='DOP853',rtol=1e-13,atol=1e-15,dense_output=True,max_step=.02)
    if not sol.success:raise RuntimeError(sol.message)
    def positions(t):
        shape=np.shape(t)
        return sol.sol(t).reshape(len(y0),6,*shape)[:,:3].swapaxes(1,-1)
    return positions

def solar_disks(pos,earth):
    x=(pos-earth)*AU
    distance=np.linalg.norm(x,axis=-1)
    to_earth=-x/distance[...,None]
    to_sun=-pos/np.linalg.norm(pos,axis=-1)[...,None]
    sep=np.arctan2(np.linalg.norm(np.cross(to_earth,to_sun),axis=-1),np.sum(to_earth*to_sun,axis=-1))
    er=np.arcsin(6371./distance)
    sr=np.arcsin(695700./(AU*np.linalg.norm(pos,axis=-1)))
    return sep,er,sr

def lost(sep,er,sr):
    a=np.clip((sep**2+er**2-sr**2)/(2*sep*er),-1,1)
    b=np.clip((sep**2+sr**2-er**2)/(2*sep*sr),-1,1)
    rad=(-sep+sr+er)*(sep+sr-er)*(sep-sr+er)*(sep+sr+er)
    return np.clip((er**2*np.arccos(a)+sr**2*np.arccos(b)-.5*np.sqrt(np.maximum(rad,0)))/(np.pi*sr**2),0,1)

def observed(times,positions,spl):
    tr=tdb(np.asarray(times))
    # Geocentric light time; topocentric correction <= Earth radius / c.
    geo=positions(tr)-spl[399](tr)
    light=np.linalg.norm(geo[0],axis=-1)*AU/299792.458/86400
    te=tr-light
    return solar_disks(positions(te),spl[399](te))

def main():
    elems,cov=orbit();spl=planets()
    z=norm.ppf(qmc.Sobol(d=6,scramble=True,seed=20241211).random_base2(8))
    samples=np.vstack([elems,elems+z@np.linalg.cholesky(cov).T])
    positions=integrate(samples,spl)
    # Compare the propagated nominal asteroid to ESA's direct geocentric astrometric output.
    angle=np.deg2rad(84381.448/3600)
    c,s=np.cos(angle),np.sin(angle)
    rotation=np.array([[1,0,0],[0,c,-s],[0,s,c]])
    errors=[]
    for line in (HERE/'esa_ephem.txt').read_text().splitlines()[9:]:
        row=line.split()
        if not row or row[0]!='11':continue
        tt=float(row[4])-60655+69.184/86400
        geo=positions(tt)[0]-spl[399](tt)
        te=tt-np.linalg.norm(geo)*AU/299792.458/86400
        direction=rotation@(positions(te)[0]-spl[399](tt))
        direction/=np.linalg.norm(direction)
        ra=np.deg2rad(15*(float(row[5])+float(row[6])/60+float(row[7])/3600))
        de=np.deg2rad(np.sign(float(row[8]))*(abs(float(row[8]))+float(row[9])/60+float(row[10])/3600))
        unit=np.array([np.cos(ra)*np.cos(de),np.sin(ra)*np.cos(de),np.sin(de)])
        errors.append(np.rad2deg(np.arctan2(np.linalg.norm(np.cross(direction,unit)),direction@unit))*3600)
    print('ESA direct sky positions, sample size/median/max arcsec',len(errors),np.median(errors),max(errors))
    def sec(h,m,s):return 3600*h+60*m+s
    def f(s,inner): 
        sep,er,sr=observed(s,positions,spl)
        return (sep-er+(sr if inner else -sr))[0]
    for bracket,inner in [((18,19),False),((18,19),True),((19,20),True),((19,20),False)]:
        k=brentq(lambda s:f(s,inner),*(x*3600 for x in bracket))
        print('UTC reception contact',utc0+timedelta(seconds=k))
    for t in [sec(18,27,41.1),sec(18,32,58.9)]:
        times=t+np.linspace(-109.5,109.5,221)
        sep,er,sr=observed(times,positions,spl)
        losses=lost(sep,er,sr)
        avg=np.trapezoid(losses,times,axis=1)/219
        print('loss nominal/clone percentiles (2.5/50/97.5)',avg[0],np.percentile(avg[1:],[2.5,50,97.5]),'center',losses[0,110])
    tt=tdb(sec(18,32,58.9));x=(positions(tt)[0]-spl[399](tt))*AU
    print('earth center distance km',np.linalg.norm(x))

if __name__=='__main__': main()
