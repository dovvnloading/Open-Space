"""Measure conditional linkage of synthetic hyperbolae at real F51 sky/time proxies."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import beta

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import link

SEED = 317801
MAX_SEPARATION = 1.5
# A larger angle is used only in the inexpensive Taylor prefilter.
PREFILTER = 2.0


def dense_path(state):
    def gravity(t,y):
        r=y[:3]
        return np.r_[y[3:],-link.MU*r/np.linalg.norm(r)**3]
    kwargs=dict(method="DOP853",rtol=1e-11,atol=1e-13,dense_output=True)
    past=solve_ivp(gravity,(0,-1.),state,**kwargs)
    future=solve_ivp(gravity,(0,8.),state,**kwargs)
    if not past.success or not future.success:
        raise ValueError("numerical trajectory failed")
    def position(t):
        out=np.empty((len(t),3))
        for mask,path in [(t<0,past),(t>=0,future)]:
            if np.any(mask):
                out[mask]=path.sol(t[mask])[:3].T
        return out
    return position


def directions(trajectory,anchor,rho,group):
    sites=np.array([link.observer(t) for t in group.t])
    dt=group.t-anchor.mid
    emission=dt.copy()
    for _ in range(3):
        emission=dt+(rho-np.linalg.norm(trajectory(emission)-sites,axis=1))/link.C
    return link.unit(trajectory(emission)-sites)


def max_sep(u,v):
    chord=np.linalg.norm(u-v,axis=1)
    return float(np.rad2deg(2*np.arcsin(min(1,chord.max()/2))))


def geometric_trial(rng, groups, index, mids, sites, sky):
    anchor=groups[index]
    rho=np.exp(rng.uniform(np.log(.45),np.log(2.8)))
    speed=rng.uniform(15.,60.)
    pos=anchor.site+rho*anchor.sky
    vel=link.unit(rng.normal(size=3))*np.sqrt((speed*link.KMS)**2+2*link.MU/np.linalg.norm(pos))
    state=np.r_[pos,vel]
    lo,hi=np.searchsorted(mids,[anchor.mid+.7,anchor.mid+6.5])
    if hi-lo<2:
        return None
    dt=mids[lo:hi]-anchor.mid
    approx=link.unit(link.short_arc(pos,vel,dt[:,None])-sites[lo:hi])
    hits=np.flatnonzero(np.linalg.norm(approx-sky[lo:hi],axis=1)<=2*np.sin(np.deg2rad(PREFILTER)/2))+lo
    pairs=[(i,j) for i in hits if mids[i]-anchor.mid<=3.25 for j in hits
           if .7<=mids[j]-mids[i]<=3.25 and groups[i].night!=groups[j].night
           and groups[i].night!=anchor.night]
    if not pairs:
        return None
    path=dense_path(state)
    witness_ids={index}|{i for pair in pairs for i in pair}
    predicted={i:directions(path,anchor,rho,groups[i]) for i in witness_ids}
    distance={i:max_sep(predicted[i],groups[i].u) for i in witness_ids}
    score,i,j=min((max(distance[index],distance[a],distance[b]),a,b) for a,b in pairs)
    if score>MAX_SEPARATION:
        return None
    ids=[index,i,j]
    return dict(anchor=anchor.name,witnesses=[groups[k].name for k in ids],indices=ids,
                max_separation_deg=score,per_night_deg=[distance[k] for k in ids],
                rho_au=float(rho),vinf_kms=float(speed),radial_helio_kms=float(pos@vel/np.linalg.norm(pos)/link.KMS),
                state=state.tolist())


def score_in_real_background(row,groups,real,pair_cache,rng,trial):
    geometry=[groups[i] for i in row["indices"]]
    path=dense_path(np.array(row["state"]))
    injected=[]
    for j,g in enumerate(geometry):
        u=directions(path,geometry[0],row["rho_au"],g)
        u=link.unit(u+rng.normal(0,.2*link.ARC,u.shape))
        tr=link.attributable(f"point-{trial}-{j}",g.t,u)
        if tr is None:
            return False,False,0
        injected.append(tr)
    pool=sorted(real+injected,key=lambda a:a.mid)
    names={a.name for a in injected}
    recovered=False
    mixed=0
    for i,a in enumerate(pool):
        if a.name in names:
            links=link.pair_matches(a,pool[i+1:])
            required=None
        else:
            if not any(.7<=z.mid-a.mid<=6.5 for z in injected):
                continue
            added=link.pair_matches(a,injected)
            if not added:
                continue
            links=pair_cache[a.name]+added
            required=names
        for trio,guess in link.triples(a,links,required):
            fit=link.correct(trio,guess)
            if fit is None:
                continue
            if {a.name for a in trio}==names:
                recovered=True
            else:
                mixed+=1
    return True,recovered,mixed


def interval(k,n):
    if not n:
        return None
    return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,
            float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,required=True)
    p.add_argument("--states",type=int,default=20000)
    p.add_argument("--fit-cap",type=int,default=100)
    args=p.parse_args()
    link.SURVEY="F51-all"
    source=args.data/"f51-all.mpc"
    real,groups,_,_=link.parse(source)
    groups=sorted(groups,key=lambda a:a.mid)
    mids=np.array([a.mid for a in groups])
    sky=np.array([a.sky for a in groups])
    sites=np.array([a.site for a in groups])
    rng=np.random.default_rng(SEED)
    geometry=[]
    for trial in range(args.states):
        row=geometric_trial(rng,groups,rng.integers(len(groups)),mids,sites,sky)
        if row is not None:
            row["trial"]=trial
            geometry.append(row)
    # Fit a reproducible uniform subsample if many states meet the fixed proxy.
    sample_rng=np.random.default_rng(SEED+1)
    which=sorted(sample_rng.choice(len(geometry),min(len(geometry),args.fit_cap),replace=False).tolist())
    real_proposed,real_pass,pairs=link.scan(real)
    for i in which:
        row=geometry[i]
        noise=np.random.default_rng(np.random.SeedSequence([SEED+2,row["trial"]]))
        eligible,recovered,mixed=score_in_real_background(row,groups,real,pairs,noise,row["trial"])
        row.update(fitted=True,eligible=eligible,recovered=recovered,mixed_passing=mixed)
    out={}
    for radius in [.1,.5,1.5]:
        g=[r for r in geometry if r["max_separation_deg"]<=radius]
        fitted=[r for r in g if r.get("fitted")]
        recovered=sum(r["recovered"] for r in fitted)
        incoming=[r for r in fitted if r["radial_helio_kms"]<0]
        out[str(radius)] = dict(geometric_states=len(g),all_anchor_states=args.states,
             geometric_fraction_ci95=interval(len(g),args.states),fitted=len(fitted),
             eligible=sum(r["eligible"] for r in fitted),recovered=recovered,
             conditional_recovery_ci95=interval(recovered,len(fitted)),
             mixed_passing=sum(r["mixed_passing"] for r in fitted),inbound=len(incoming),
             inbound_recovered=sum(r["recovered"] for r in incoming))
    summary=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),seed=SEED,
                 fit_cap=args.fit_cap,usable_anchors=len(groups),utc_nights=len({a.night for a in groups}),
                 anchor_state_trials=args.states,real_proposed=real_proposed,real_passing=real_pass,by_radius_deg=out)
    # All geometrically supported states are retained, including those not sampled for fitting.
    for row in geometry:
        del row["indices"]
    (args.data/"pointing-trials.json").write_text(json.dumps(geometry,indent=2)+"\n")
    (args.data/"pointing-result.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("POINTING_RESULT="+json.dumps(summary))

if __name__=="__main__":
    main()
