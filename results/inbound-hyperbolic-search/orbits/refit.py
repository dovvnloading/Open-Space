"""Independent sky-plane fits, with the barycentric inbound branch followed in DE440s.

Run acquire.py first. Outputs are conditional on the retained ground optical
records and each force law. Statistics here are sensitivity diagnostics.
"""
from collections import Counter
from pathlib import Path
import base64
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import re
import sys
import traceback
import urllib.parse

import erfa
from jplephem.spk import SPK
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

from acquire import CACHE, DATA, NAMES, retrieve

# Same eight planetary systems, solar GM, frame orientation conventions and kernel
# as the repository's earlier independent integration. Its single-object state
# initializer is specific to the 704021 encounter.
prior_path = Path(__file__).resolve().parents[2] / "borisov-704021" / "propagate.py"
spec = importlib.util.spec_from_file_location("earlier_propagation", prior_path)
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
AU, DAY = prior.AU, prior.DAY
MUS = np.r_[prior.GM, prior.PLANET_GM]
GM_TOTAL = MUS.sum()
ARCSEC = 180 * 3600 / np.pi
C = 173.144632674240
OBS_TIME_LIMIT = 7  # Light-time interpolation pad before first reception in days.

def stations():
    raw = (CACHE / "stations.txt").read_text(errors="replace")
    rows = re.sub(r"<[^>]*>","",raw).splitlines()
    out = {}
    for row in rows:
        try:
            lam, rho_c, rho_s = float(row[4:13]), float(row[13:21]), float(row[21:30])
        except ValueError:
            continue
        lam = np.deg2rad(lam)
        out[row[:3]] = (6378.137 / AU) * np.array([rho_c*np.cos(lam),rho_c*np.sin(lam),rho_s])
    return out

def observations(obj, eph):
    sites = stations()
    groups, jd, angles, geoc, seen = [], [], [], [], set()
    rejected = Counter()
    lines = (DATA/(obj+".obs")).read_text().splitlines()
    for row in lines:
        mode, site = row[14], row[77:80]
        if mode not in " CBAOPTMc" or site not in sites:
            rejected["unsupported mode or station"] += 1
            continue
        try:
            year, month, day = row[15:32].split()
            y,m,d = int(year),int(month),float(day)
            t = sum(erfa.cal2jd(y,m,1)) + d-1
            h,minute,s = map(float,row[32:44].split())
            deg,amin,sec = map(float,row[45:56].split())
            ra = 15*(h+minute/60+s/3600)
            dec = (-1 if row[44]=="-" else 1)*(deg+amin/60+sec/3600)
        except (ValueError, IndexError):
            rejected["unreadable"]+=1
            continue
        key=(site,t,ra,dec)
        if key in seen:
            rejected["duplicate"]+=1
            continue
        seen.add(key)
        groups.append((site,int(t-.5)))
        jd.append(t)
        angles.append([ra,dec])
        geoc.append(sites[site])
    if len(jd)<7:
        raise ValueError(f"Only {len(jd)} supported optical records")
    utc=np.asarray(jd)
    tai=erfa.utctai(2400000.5,utc-2400000.5)
    tt=erfa.taitt(*tai)
    # The 1877 catalogue dates approximate UT1. UTC leap seconds did not exist.
    # Espenak & Meeus 1860–1900 delta-T polynomial (year since 1860).
    historical=(utc>=2400410.5) & (utc<2415020.5)
    if np.any(historical):
        yy,mm,dd,frac=erfa.jd2cal(2400000.5,utc[historical]-2400000.5)
        y=yy+(mm-.5)/12-1860
        delta=7.62+.5737*y-.251754*y**2+.01680668*y**3-.0004473624*y**4+y**5/233174
        tt=(np.full_like(utc,2400000.5),np.where(historical,utc-2400000.5,tt[0]+tt[1]-2400000.5))
        tt[1][historical]+=delta/DAY
    tdb=tt[0]+tt[1]  # TT-TDB <2 ms
    rot=erfa.c2t06a(*tt,np.full_like(utc,2400000.5),utc-2400000.5,0.,0.)
    geoc=np.einsum("nji,nj->ni",rot,geoc)
    earth=(eph[0,3].compute(tdb)+eph[3,399].compute(tdb)).T/AU
    # Post-1960 weighting uses one arcsec per sky coordinate; earlier optical
    # measurements get a five arcsec floor. Limit station-night total weight.
    counts=Counter(groups)
    sigma=np.array([np.sqrt(max(1,counts[g]/4)) for g in groups])*np.where(utc<2436934.5,5.,1.)
    ra,dec=np.deg2rad(angles).T
    p=np.column_stack([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
    e=np.column_stack([-np.sin(ra),np.cos(ra),np.zeros_like(ra)])
    n=np.cross(p,e)
    return tdb,earth+geoc,e,n,sigma,groups,dict(rejected)

def seed(obj, epoch, receipts):
    name=obj+";" if obj.endswith("I") else "DES="+obj.replace("-"," ")+";"
    params={"format":"json","COMMAND":"'"+name+"'","OBJ_DATA":"NO","EPHEM_TYPE":"VECTORS",
            "CENTER":"'500@0'","TLIST":f"'{epoch:.8f}'","CSV_FORMAT":"YES","REF_PLANE":"FRAME",
            "VEC_TABLE":"2","OUT_UNITS":"AU-D","VEC_CORR":"NONE"}
    url="https://ssd.jpl.nasa.gov/api/horizons.api?"+urllib.parse.urlencode(params)
    raw=retrieve(url,obj+"-seed.json",receipts)
    doc=json.loads(raw)
    if "$$SOE" not in doc.get("result",""):
        raise ValueError("Horizons seed rejected: "+str(doc)[:240])
    row=next(csv.reader(io.StringIO(doc["result"].split("$$SOE")[1].split("$$EOE")[0].strip())))
    return np.array([float(x) for x in row[2:8]])

def ephemeris(eph, first,last,step):
    days=np.arange(np.floor(first)-1,np.ceil(last)+step+1,step)
    xyz=np.stack([eph[0,j].compute(days).T/AU for j in [10,1,2,3,4,5,6,7,8]],axis=1)
    return CubicSpline(days,xyz,axis=0,extrapolate=False)

def derivative(t,y,positions,npar):
    x=y[:3]
    d=positions(t)-x
    r2=np.einsum("ij,ij->i",d,d)
    inv=r2**-1.5
    mus=MUS.copy()
    if npar==7:
        # A1 in units 10^-8 au/day^2, positive outward, inverse square.
        mus[0]-=y[6]*1e-8
    coeff=mus*inv
    accel=coeff@d
    vmat=y[npar:].reshape(6,npar)
    grad=3*np.einsum("i,ij,ik->jk",coeff/r2,d,d)-coeff.sum()*np.eye(3)
    partial=grad@vmat[:3]
    if npar==7:
        partial[:,-1]-=1e-8*inv[0]*d[0]
    if npar==8:
        # Marsden, Sekanina & Yeomans 1973, symmetric water-ice law normalized at 1 au.
        r=np.sqrt(r2[0])
        er=-d[0]/r
        vhel=y[3:6]-positions(t,1)[0]
        vr=er@vhel
        vt=vhel-vr*er
        speed=np.linalg.norm(vt)
        et=vt/speed
        z=(r/2.808)**5.093
        g=.111262*(r/2.808)**(-2.15)*(1+z)**(-4.6142)
        dg_dr=g*(-2.15-4.6142*5.093*z/(1+z))/r
        per=np.eye(3)-np.outer(er,er)
        det=(np.eye(3)-np.outer(et,et))/speed
        der=per/r
        det_dx=det@(-np.outer(er,vhel@der)-vr*der)
        a1,a2=y[6:8]*1e-8
        direction=a1*er+a2*et
        accel+=g*direction
        px=g*(a1*der+a2*det_dx)+np.outer(direction,dg_dr*er)
        pv=g*a2*(det@per)
        partial+=px@vmat[:3]+pv@vmat[3:]
        partial[:,6:8]+=1e-8*g*np.column_stack([er,et])
    base=np.r_[y[3:6],accel,np.zeros(npar-6)]
    return np.r_[base,vmat[3:].ravel(),partial.ravel()]

def integrate(params,start,stop,positions,events=None,dense=True):
    n=len(params)
    aug=np.r_[params,np.eye(6,n).ravel()]
    solution=solve_ivp(derivative,(start,stop),aug,method="DOP853",rtol=1e-11,atol=2e-13,
                       args=(positions,n),dense_output=dense,events=events)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution

def fit(seed_state,epoch,t,observers,east,north,sig,positions,outward=False,return_predictor=False):
    npar=len(seed_state)
    memo={}
    if len(t)<max(7,npar+1):
        raise ValueError("Too few positions for independent refit")
    def calc(state):
        key=tuple(state)
        if key not in memo:
            memo.clear()
            sol=integrate(state,epoch,t.min()-OBS_TIME_LIMIT,positions).sol
            emit=t.copy()
            for _ in range(4):
                y=sol(emit)
                delta=y[:3].T-observers
                distance=np.linalg.norm(delta,axis=1)
                emit=t-distance/C
            y=sol(emit)
            d=y[:3].T-observers
            dist=np.linalg.norm(d,axis=1)
            u=d/dist[:,None]
            vel=y[3:6].T
            deriv=y[npar:].T.reshape(-1,6,npar)[:,:3]
            tau=np.einsum("ni,nij->nj",u,deriv)/(C+np.sum(u*vel,axis=1))[:,None]
            deriv=deriv-vel[:,:,None]*tau[:,None,:]
            du=(deriv-u[:,:,None]*np.einsum("ni,nij->nj",u,deriv)[:,None,:])/dist[:,None,None]
            directions=np.stack([east,north],axis=1)
            resid=ARCSEC*np.einsum("nij,nj->ni",directions,u)
            jmat=ARCSEC*np.einsum("nik,nkj->nij",directions,du)
            memo[key]=resid,jmat
        return memo[key]
    mask=np.ones(len(t),dtype=bool)
    def f(x):
        return (calc(x)[0][mask]/sig[mask,None]).ravel()
    def j(x):
        return (calc(x)[1][mask]/sig[mask,None,None]).reshape(-1,npar)
    lower=np.full(npar,-np.inf)
    if outward:
        lower[6]=0
    x0=seed_state
    for passno in range(3):
        result=least_squares(f,x0,jac=j,x_scale="jac",bounds=(lower,np.inf),loss="soft_l1" if passno==0 else "linear",
                             ftol=1e-8,xtol=1e-10,gtol=1e-6,max_nfev=70)
        if not result.success: raise ValueError(result.message)
        x0=result.x
        if passno<2:
            mask=np.linalg.norm(calc(x0)[0],axis=1)<4*sig
            if mask.sum()<max(7,.5*len(mask)):
                raise ValueError("Fit excludes too many observations")
    res=calc(x0)[0]
    u,s,vt=np.linalg.svd(j(x0),full_matrices=False)
    if s[-1]<=s[0]*1e-14:
        raise ValueError("Singular normal matrix")
    redchi=np.sum(f(x0)**2)/(2*mask.sum()-npar)
    cov=(vt.T/s**2)@vt*max(redchi,1)
    result=(x0,cov,mask,res,redchi)
    return (*result,calc) if return_predictor else result

def inbound(obj,params,cov,epoch,eph):
    # The outer-boundary state is measured with respect to SSB, using total
    # system GM. C/2024 L5 uses a deliberate 2020 pre-Saturn encounter epoch.
    stop = 2458849.5 if obj=="2024-L5" else max(max(seg.start_jd for seg in eph.segments if seg.center==0)+12,epoch-450*365.25)
    positions=ephemeris(eph,stop-2,epoch+1,1.)
    def crossing(radius,terminal):
        def event(t,y,*args):
            return np.linalg.norm(y[:3])-radius
        event.direction=1
        event.terminal=terminal
        return event
    events=[crossing(50,False),crossing(250,True)]
    sol=integrate(params,epoch,stop,positions,events,dense=False)
    summaries={}
    n=len(params)
    for label,when,state in [("comparison_epoch",sol.t[-1],sol.y[:,-1])]+[
         (str(radius)+"au",ts[0],ys[0]) for radius,ts,ys in zip([50,250],sol.t_events,sol.y_events) if len(ts)]:
        pos,vel=state[:3],state[3:6]
        radius=np.linalg.norm(pos)
        ir=2/radius-vel@vel/GM_TOTAL
        g=np.r_[-2*pos/radius**3,-2*vel/GM_TOTAL]
        stm=state[n:].reshape(6,n)
        vr=pos@vel/radius
        if label.endswith("au"):
            # Different solutions reach this sphere at different times.
            tgrad=-(pos/radius)@stm[:3]/vr
            stm=stm+np.outer(derivative(when,state,positions,n)[:6],tgrad)
        sig=np.sqrt(max(0,(g@stm)@cov@(g@stm)))
        summaries[label]={"jd_tdb":when,"radius_au":radius,"inbound_radial_au_d":vr,
                          "inverse_a_au_inv":ir,"formal_sigma_inverse_a":sig,
                          "v_inf_km_s":np.sqrt(-ir*GM_TOTAL)*AU/DAY if ir<0 else None}
    return summaries

def parabolic_profile(epoch,eph,nominal,laws,mask,sigma,evaluate):
    """Fit the same retained positions while imposing zero energy at inbound 50 au.

    A tight energy residual is added only to find the constrained solution;
    the reported delta chi-square uses only the actual sky measurements.
    """
    pos=ephemeris(eph,epoch-130*365.25,epoch+1,1.)
    n=len(nominal)
    cache={}
    def event(t,y,*args):
        return np.linalg.norm(y[:3])-50
    event.direction=1
    event.terminal=True
    def alpha(x):
        key=tuple(x)
        if key not in cache:
            cache.clear()
            sol=integrate(x,epoch,epoch-129*365.25,pos,[event],False)
            if not len(sol.t_events[0]):
                raise ValueError("Trial missed inbound 50 au within 129 years")
            when=sol.t_events[0][0]; y=sol.y_events[0][0]
            r,v=y[:3],y[3:6]; R=np.linalg.norm(r); vr=r@v/R
            if vr>=0:
                raise ValueError("Wrong branch at 50 au")
            stm=y[n:].reshape(6,n)
            tgrad=-(r/R)@stm[:3]/vr
            stm=stm+np.outer(derivative(when,y,pos,n)[:6],tgrad)
            g=np.r_[-2*r/R**3,-2*v/GM_TOTAL]@stm
            cache[key]=(2/R-v@v/GM_TOTAL,g)
        return cache[key]
    def measurements(x):
        res,jac=evaluate(x)
        return (res[mask]/sigma[mask,None]).ravel(),(jac[mask]/sigma[mask,None,None]).reshape(-1,n)
    tol=2e-5 # au^-1. The achieved constraint error is reported.
    def f(x): return np.r_[measurements(x)[0],alpha(x)[0]/tol]
    def j(x): return np.vstack([measurements(x)[1],alpha(x)[1]/tol])
    lo=np.full(n,-np.inf)
    if laws=="marsden_outward": lo[6]=0
    out=least_squares(f,nominal,jac=j,bounds=(lo,np.inf),x_scale="jac",max_nfev=120,
                      ftol=1e-9,xtol=1e-11,gtol=1e-6)
    if not out.success: raise ValueError(out.message)
    before=measurements(nominal)[0]; after=measurements(out.x)[0]
    return {"target_inverse_a_at_inbound_50au":0.,"achieved":alpha(out.x)[0],
            "delta_chi2_sky_only":float(after@after-before@before),"same_positions":int(mask.sum()),
            "profile_state_and_accelerations":out.x.tolist(),"n_function_evals":out.nfev}

def check_holdouts(obj,initial,epoch,eph,t,observers,east,north,sig,groups,positions):
    """Drop endpoint nights and the top-contributing observatory; evaluate withheld sky positions."""
    groups=list(groups)
    nights=np.array([night for site,night in groups])
    station=np.array([site for site,night in groups])
    top=Counter(station).most_common(1)[0][0]
    cuts={"first_utc_night":nights==nights.min(),"last_utc_night":nights==nights.max(),
          "largest_station_"+str(top):station==top}
    if np.any(station=="I41"):
        cuts["ZTF_I41"]=station=="I41"
    tests={}
    for name,drop in cuts.items():
        try:
            x,cov,mask,res,redchi=fit(initial,epoch,t[~drop],observers[~drop],east[~drop],north[~drop],
                                      sig[~drop],positions)
            path=integrate(x,epoch,t.min()-OBS_TIME_LIMIT,positions).sol
            emit=t[drop].copy()
            for _ in range(5):
                delta=path(emit)[:3].T-observers[drop]
                dist=np.linalg.norm(delta,axis=1)
                emit=t[drop]-dist/C
            u=delta/dist[:,None]
            error=ARCSEC*np.column_stack([np.sum(u*east[drop],axis=1),np.sum(u*north[drop],axis=1)])
            tests[name]={"withheld":int(drop.sum()),"fitted":int(mask.sum()),
                         "withheld_rms_per_coordinate_as":float(np.sqrt(np.mean(error**2))),
                         "boundaries":inbound(obj,x,cov,epoch,eph)}
        except Exception as exc:
            tests[name]={"failed":str(exc)}
    return tests

def run():
    receipts=[]
    report={}
    # Carry the exact input bytes through the read-only runner for durable archiving.
    for src in sorted(DATA.glob("*.obs")):
        packed=gzip.compress(src.read_bytes(),mtime=0)
        encoded=base64.b64encode(packed).decode()
        for i in range(0,len(encoded),6000):
            print("OBS_CHUNK",src.name,i,encoded[i:i+6000],flush=True)
    output=DATA/"refits.json"
    kernel=CACHE/"de440.bsp"
    with SPK.open(str(kernel if kernel.exists() else CACHE/"de440s.bsp")) as eph:
        for obj in NAMES:
            entry={"models":{}}
            report[obj]=entry
            try:
                t,stations_,east,north,sig,groups,dropped=observations(obj,eph)
                epoch=float(np.ceil(t.max())+1)
                initial=seed(obj,epoch,receipts)
                entry["initial_horizons_seed_ssb_icrf"]=initial.tolist()
                pos=ephemeris(eph,t.min()-OBS_TIME_LIMIT-2,epoch+1,.25)
                entry.update(n_ground_optical=len(t),n_station_nights=len(set(groups)),
                             observed_jd_tdb=[float(t.min()),float(t.max())],excluded_rows=dropped)
                models=["gravity"]
                if obj in {"2025-V1","1I","1997-P2","2008-J4"}:
                    models.append("radial_r2")
                if obj=="2025-V1":
                    models.extend(["marsden_signed","marsden_outward"])
                if obj=="1877-G2":
                    models.extend(["gravity_2as","gravity_10as"])
                gravity_seed=initial
                for kind in models:
                    try:
                        scale={"gravity_2as":.4,"gravity_10as":2.}.get(kind,1.)
                        weights=sig*scale
                        nng=2 if kind.startswith("marsden") else int(kind=="radial_r2")
                        state=np.r_[gravity_seed,np.zeros(nng)]
                        x,cov,mask,res,redchi,predict=fit(
                            state,epoch,t,stations_,east,north,weights,pos,
                            outward=kind=="marsden_outward",return_predictor=True)
                        if kind=="gravity":
                            gravity_seed=x.copy()
                        summary={"epoch_jd_tdb":epoch,"state_ssb_icrf":x.tolist(),"covariance":cov.tolist(),
                                 "used":int(mask.sum()),"used_indices_in_supported_records":np.flatnonzero(mask).tolist() if obj in ("1877-G2","2008-J4") else None,
                                 "omitted_by_clipping":int((~mask).sum()),
                                 "rms_per_coordinate_as":float(np.sqrt(np.mean(res[mask]**2))),
                                 "weighted_reduced_chi2":float(redchi),
                                 "boundaries":inbound(obj,x,cov,epoch,eph)}
                        if obj=="2025-V1":
                            try:
                                summary["zero_energy_profile"]=parabolic_profile(epoch,eph,x,kind,mask,sig,predict)
                            except Exception as exc:
                                summary["zero_energy_profile"]={"failed":str(exc)}
                        if obj=="2025-V1" and kind in {"gravity","marsden_signed"}:
                            summary["holdouts"]=check_holdouts(obj,x,epoch,eph,t,stations_,east,north,sig,groups,pos)
                        entry["models"][kind]=summary
                        print("FIT",obj,kind,summary["used"],summary["rms_per_coordinate_as"],
                              json.dumps(summary["boundaries"]),flush=True)
                    except Exception as exc:
                        entry["models"][kind]={"failed":str(exc)}
                        traceback.print_exc()
            except Exception as exc:
                entry["failed"]=str(exc)
                traceback.print_exc()
            output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
            print("OBJECT_OUT",obj,json.dumps(entry),flush=True)
            (DATA/"initialization_queries.json").write_text(json.dumps(receipts,indent=2)+"\n")
    print("SEED_QUERIES",json.dumps(receipts),flush=True)
    if not any("boundaries" in m for entry in report.values() for m in entry["models"].values()):
        raise RuntimeError("No independent orbit propagated")

if __name__=="__main__":
    run()
