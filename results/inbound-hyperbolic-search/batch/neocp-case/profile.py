"""Independent short-arc range profile for CEXJ3W2. No catalogue orbit is used."""
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import erfa
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
AU = 149597870.7
DAY = 86400.
KMS = DAY/AU
ARC = np.deg2rad(1/3600)
C = 299792.458*KMS
GM_S = 0.0002959122082855911
GM_E = 398600.435507*DAY**2/AU**3
GM_M = 4902.800118*DAY**2/AU**3
SITES = {"G96": (249.21128,0.845107,0.533611),
         "I52": (249.21108,0.845109,0.533609)}
JD0 = 2461342.5
TT_UTC = 69.184/DAY
IERS = "https://maia.usno.navy.mil/ser7/finals2000A.all"


def unit(a):
    return a/np.linalg.norm(a,axis=-1,keepdims=True)


def gravity(x):
    return x/np.linalg.norm(x,axis=-1,keepdims=True)**3


def eop():
    cache=HERE/"evidence/eop.json"
    if cache.is_file():
        return json.loads(cache.read_text())
    with urlopen(IERS,timeout=40) as r:
        raw=r.read()
        retrieved=r.headers.get("Date")
    lines=[l for l in raw.decode("ascii").splitlines() if l[7:15].strip() in {"61342.00","61343.00"}]
    if len(lines)!=2:
        raise ValueError("Daily Earth orientation measurements unavailable")
    return {"url":IERS,"date":retrieved,"full_source_sha256":hashlib.sha256(raw).hexdigest(),
            "raw_rows":lines,"values":[[float(l[18:27]),float(l[37:46]),float(l[58:68])] for l in lines]}


def site(station,frac,eops):
    lon,rc,rs=SITES[station]
    lon=np.deg2rad(lon)
    itrf=6378.137/AU*np.array([rc*np.cos(lon),rc*np.sin(lon),rs])
    xp,yp,ut1=eops[0]+frac*(eops[1]-eops[0])
    c2t=erfa.c2t06a(JD0,frac+TT_UTC,JD0,frac+ut1/DAY,xp*ARC,yp*ARC)
    return c2t.T@itrf


class Model:
    def __init__(self,lines,eops):
        self.t=np.array([float(l[23:32])-29 for l in lines])
        self.epoch=float(np.mean(self.t))
        self.t-=self.epoch
        self.stations=[l[77:80] for l in lines]
        self.flags=np.array([l[64]=="x" for l in lines])
        sky=[]
        for l in lines:
            h,m,s=map(float,l[32:44].split())
            d,a,b=map(float,l[45:56].split())
            ra=np.deg2rad(15*(h+m/60+s/3600))
            dec=np.deg2rad((1 if l[44]=="+" else -1)*(d+a/60+b/3600))
            sky.append([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
        self.sky=np.array(sky)
        self.u=unit(self.sky.mean(axis=0))
        self.e=unit(np.cross([0,0,1],self.u))
        self.n=np.cross(self.u,self.e)
        rates=self.t@self.sky/(self.t@self.t)
        self.rates=np.array([rates@self.e,rates@self.n])
        self.site0=site("G96",self.epoch,eops)
        self.sitev=(site("G96",self.epoch+.0001,eops)-site("G96",self.epoch-.0001,eops))/.0002
        self.sites=np.array([site(o,t+self.epoch,eops) for o,t in zip(self.stations,self.t)])
        ph,pb=erfa.epv00(JD0,self.epoch+TT_UTC)
        self.heliop,self.heliov=ph["p"],ph["v"]
        self.baryv=pb["v"]
        moon=erfa.moon98(JD0,self.epoch+TT_UTC)
        self.moonp,self.moonv=moon["p"],moon["v"]
        self.eartha=-GM_S*gravity(self.heliop)+GM_M*gravity(self.moonp)
        self.baseline_km=float(np.linalg.norm(site("G96",self.epoch,eops)-site("I52",self.epoch,eops))*AU)

    def state(self,x,rho):
        # x contains tangent offsets (arcsec), angular-rate corrections (arcsec/day),
        # and topocentric range rate (km/s) at a common reception-time epoch.
        u=unit(self.u+ARC*(x[0]*self.e+x[1]*self.n))
        tang=self.e*(self.rates[0]+x[2]*ARC)+self.n*(self.rates[1]+x[3]*ARC)
        return self.site0+rho*u,self.sitev+x[4]*KMS*u+rho*(tang-(tang@u)*u)

    def acc(self,t,p):
        sun=-self.heliop-t[...,None]*self.heliov
        moon=self.moonp+t[...,None]*self.moonv
        return -GM_E*gravity(p)+GM_S*(gravity(sun-p)-gravity(sun))+GM_M*(gravity(moon-p)-gravity(moon))

    def propagate(self,pos,vel,t,nstep):
        # Integrate geocentric inertial trajectories simultaneously at all requested dates.
        h=t/nstep
        step=h[:,None]
        x=np.broadcast_to(pos,(len(t),3)).copy()
        v=np.broadcast_to(vel,(len(t),3)).copy()
        for k in range(nstep):
            a1=self.acc(k*h,x)
            a2=self.acc((k+.5)*h,x+step*v/2)
            a3=self.acc((k+.5)*h,x+step*(v+step*a1/2)/2)
            a4=self.acc((k+1)*h,x+step*(v+step*a2/2))
            x+=step*(v+step*(a1+a2+a3)/6)
            v+=step*(a1+2*a2+2*a3+a4)/6
        return x

    def residual(self,x,rho,nstep=2):
        pos,vel=self.state(x,rho)
        # Translate from geocentric emission coordinates to the reception-time observer.
        te=self.t-np.linalg.norm(pos+vel*self.t[:,None]-self.sites,axis=1)/C
        for _ in range(2):
            p=self.propagate(pos,vel,te,nstep)
            delay=(te-self.t)[:,None]
            earth_shift=delay*(self.baryv+((te+self.t)/2)[:,None]*self.eartha)
            delta=p+earth_shift-self.sites
            te=self.t-np.linalg.norm(delta,axis=1)/C
        dif=unit(delta)-self.sky
        return np.column_stack((dif@self.e,dif@self.n))/ARC

    def elements(self,x,rho):
        p,v=self.state(x,rho)
        r_s,v_s=self.heliop+p,self.heliov+v
        h=np.cross(r_s,v_s)
        ecc=np.linalg.norm(np.cross(v_s,h)/GM_S-unit(r_s))
        energy=(v_s@v_s-2*GM_S/np.linalg.norm(r_s))
        return dict(e=float(ecc),inv_a_au=float(-energy/GM_S),specific_v2_excess_kms2=float(energy/KMS**2),
                    helio_vinf_kms=float(np.sqrt(energy)/KMS) if energy>0 else None,
                    geocentric_km=float(np.linalg.norm(p)*AU),topo_rate_kms=float(x[4]),
                    geocentric_position_au=p.tolist(),geocentric_velocity_au_day=v.tolist())


def profile(mod,mask,distances):
    previous=np.zeros(5)
    rows=[]
    for rho in distances:
        # Exclude line-of-sight starting guesses that cross the observer in this arc.
        rate_limit=min(400.,.85*rho/(np.max(np.abs(mod.t))*KMS))
        lo=np.array([-700.,-700.,-40000.,-40000.,-rate_limit])
        last=np.clip(previous,lo+.0001,-lo-.0001)
        def resid(x):
            return mod.residual(x,rho)[mask].ravel()
        opt=least_squares(resid,last,bounds=(lo,-lo),x_scale=[5,5,100,100,10],
                          xtol=5e-9,ftol=1e-9,gtol=5e-7,max_nfev=65)
        previous=opt.x
        # Finer integration is an independent numerical-convergence check for each row.
        fine=mod.residual(opt.x,rho,8)
        co=mod.residual(opt.x,rho,2)
        score=float(np.sum(fine[mask]**2))
        rows.append({**mod.elements(opt.x,rho),"rho_au":float(rho),"rho_km":float(rho*AU),
                     "rms_radial_arcsec":float(np.sqrt(score/np.sum(mask))),"sum_squares_arcsec2":score,
                     "profile_n":int(np.sum(mask)),"converged":bool(opt.success),
                     "numerical_max_delta_arcsec":float(np.max(np.linalg.norm(fine-co,axis=1))),
                     "parameters":opt.x.tolist(),"all_residuals_arcsec":fine.tolist()})
    return rows


def main():
    raw=gzip.decompress((HERE.parent/"evidence/neocp_history_initial.json.gz").read_bytes())
    source=next(t for t in json.loads(raw)["results"] if t["id"]=="CEXJ3W2")
    lines=source["observations"]
    digest=hashlib.sha256(("\n".join(lines)+"\n").encode()).hexdigest()
    assert digest=="648d3bf459e79e5bcdbea8eadc7e6f31e1b6ef2b0739064a99e0f1b2b894e796" and len(lines)==19
    earth=eop()
    mod=Model(lines,np.array(earth["values"]))
    grid=np.geomspace(.0005,.10,110)
    result={"input_sha256":digest,"input_receipt":source["receipt"],"eop":earth,
            "epoch_jd_utc":JD0+mod.epoch,"arc_minutes":float(np.ptp(mod.t)*1440),
            "station_counts":dict(collections.Counter(mod.stations)),"site_separation_km":mod.baseline_km,
            "range_grid_au":list(map(float,[grid[0],grid[-1]])),"grid_count":len(grid),"fits":{}}
    detailed=[]
    for name,mask in [("all19",np.ones(19,dtype=bool)),("published15",~mod.flags)]:
        rows=profile(mod,mask,grid)
        eligible=[r for r in rows if r["converged"] and r["numerical_max_delta_arcsec"]<.02]
        if not eligible:
            raise ValueError("No numerically reliable profile fits")
        best=min(eligible,key=lambda r:r["sum_squares_arcsec2"])
        fitted={}
        for regime,sel in [("best",eligible),("bound",[r for r in eligible if r["e"]<1]),
                           ("securely_bound_nominal",[r for r in eligible if r["specific_v2_excess_kms2"] < -100]),
                           ("strongly_unbound",[r for r in eligible if r["specific_v2_excess_kms2"]>100]),
                           ("near_parabolic",[min(eligible,key=lambda r:abs(r["e"]-1))])]:
            if sel:
                win=min(sel,key=lambda r:r["sum_squares_arcsec2"])
                fitted[regime]=win
        for allowance in [1.,4.,9.]:
            cand=[r for r in eligible if r["sum_squares_arcsec2"]<=best["sum_squares_arcsec2"]+allowance]
            fitted[f"delta_sse_{allowance:g}"]={"rho_au":[min(r["rho_au"] for r in cand),max(r["rho_au"] for r in cand)],
                    "eccentricity":[min(r["e"] for r in cand),max(r["e"] for r in cand)],"grid_points":len(cand)}
        result["fits"][name]={"n":int(np.sum(mask)),"reliable_grid_points":len(eligible),"solutions":fitted}
        detailed.extend(dict(sample=name,**{k:v for k,v in r.items() if k not in {"parameters","all_residuals_arcsec","geocentric_position_au","geocentric_velocity_au_day"}}) for r in rows)
    output=HERE/"evidence"
    output.mkdir(exist_ok=True)
    (output/"profile.json").write_text(json.dumps(result,indent=2)+"\n")
    (output/"eop.json").write_text(json.dumps(earth,indent=2)+"\n")
    with (output/"ranges.csv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=detailed[0].keys(),lineterminator="\n")
        w.writeheader()
        w.writerows(detailed)
    for name,d in result["fits"].items():
        for label in ["bound","strongly_unbound","best"]:
            r=d["solutions"].get(label,{})
            print(name,label,r.get("rho_au"),r.get("e"),r.get("rms_radial_arcsec"),flush=True)


if __name__=="__main__":
    main()
