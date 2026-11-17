"""Three-night, heliocentric trial-state linkage in a Rubin MPC1992 slice."""
import argparse
import collections
import datetime as dt
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from itertools import combinations
from pathlib import Path

import erfa
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.stats import beta

KM = 149597870.7
KMS = 86400 / KM
MU = 0.0002959122082855911
C = 299792.458 * KMS
ARC = np.pi / (180 * 3600)
RNG_SEED = 702613
SURVEY = "X05"
SITES = {"X05": (289.25058, .864981, -.500958), "F51": (203.74409, .936241, .351543)}
J2000 = dt.date(2000, 1, 1)
RANGES, RADIAL = [a.ravel() for a in np.meshgrid(np.geomspace(.4, 3, 13), np.arange(-90, 91, 15) * KMS)]


def norm(x):
    return np.linalg.norm(x, axis=-1, keepdims=True)


def unit(x):
    return x / norm(x)


@lru_cache(maxsize=200000)
def observer(t):
    # ERFA epv00 = heliocentric BCRS equatorial axes. UTC approximates UT1.
    # TT approximates TDB; polar motion is omitted. MPC X05 geocentric factors.
    tt = t + 69.184 / 86400
    pv, _ = erfa.epv00(2451545., tt)
    lon, rc, rs = SITES[SURVEY]
    lon = np.deg2rad(lon)
    xyz = np.array([rc*np.cos(lon), rc*np.sin(lon), rs]) * 6378.137/KM
    return pv["p"] + erfa.c2t06a(2451545., tt, 2451545., t, 0., 0.).T @ xyz


@dataclass
class Track:
    name: str
    t: np.ndarray
    u: np.ndarray
    mid: float
    sky: np.ndarray
    rate: np.ndarray
    night: int
    site: np.ndarray


def attributable(name, t, u):
    t, indices = np.unique(t, return_index=True)
    u = np.asarray(u)[indices]
    span = np.ptp(t)
    if len(t) < 2 or not 5/1440 <= span <= .25:
        return None
    mid = np.mean(t)
    d = t - mid
    sky = unit(u.mean(axis=0))
    rate = d @ u / (d @ d)
    rate -= (rate @ sky) * sky
    linear_error = norm(u - unit(sky + d[:, None]*rate))[:, 0] / ARC
    if linear_error.max() > 2 or not .08 <= np.rad2deg(np.linalg.norm(rate)) <= 5:
        return None
    return Track(name, t, u, mid, sky, rate, int(np.floor(mid+.5)), observer(mid))


def parse(path):
    raw = collections.defaultdict(list)
    qc = collections.Counter()
    for l in path.read_text().splitlines():
        qc["lines"] += 1
        name = l[:12].strip()
        if not name or l[14] not in "CB": # optical CCD / CMOS only
            qc["blank_or_nonoptical"] += 1
            continue
        y, m, dd = l[15:32].split()
        day = dt.date(int(y), int(m), int(float(dd)))
        t = (day-J2000).days - .5 + float(dd) % 1
        h, m, s = map(float, l[32:44].split())
        ra = np.deg2rad(15*(h + m/60 + s/3600))
        d, m, s = map(float, l[45:56].split())
        dec = np.deg2rad((1 if l[44] == "+" else -1)*(d+m/60+s/3600))
        u = [np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)]
        raw[name+"@"+str(day)].append((t, u))
    tracks = []
    for name, lines in raw.items():
        t, u = zip(*lines)
        track = attributable(name, np.array(t), u)
        if track is not None:
            tracks.append(track)
    qc["submitted_name_utc_days"] = len(raw)
    qc["eligible_tracks"] = len(tracks)
    return sorted(tracks, key=lambda a: a.mid), dict(qc)


def motion(state, times):
    # Independent numerical, solar two-body trajectory for injection and final fitting.
    t = np.asarray(times)
    def rhs(_, y):
        return np.r_[y[3:], -MU*y[:3] / np.linalg.norm(y[:3])**3]
    out = np.empty((len(t),3))
    for mask in [t < 0, t >= 0]:
        if np.any(mask):
            edge = t[mask][np.argmax(np.abs(t[mask]))]
            if abs(edge) < 1e-12:
                out[mask] = state[:3]
                continue
            f = solve_ivp(rhs, (0, edge), state, dense_output=True, method="DOP853",rtol=1e-10,atol=1e-12)
            if not f.success:
                raise ValueError("propagation failed")
            out[mask] = f.sol(t[mask])[:3].T
    return out


def short_arc(pos, vel, days):
    # Third-order heliocentric Taylor prediction only for proposing a linkage.
    r = norm(pos)
    acc = -MU*pos/r**3
    jerk = -MU*(vel/r**3 - 3*pos*np.sum(pos*vel,axis=-1,keepdims=True)/r**5)
    return pos + days*(vel + days*(acc/2 + days*jerk/6))


def initial(seed):
    pos = seed.site + RANGES[:,None]*seed.sky
    site_v = (observer(seed.mid+.0001)-observer(seed.mid-.0001))/.0002
    vel = site_v + RANGES[:,None]*seed.rate + RADIAL[:,None]*seed.sky
    return pos, vel


def proposals(seed, later):
    pos, vel = initial(seed)
    links = []
    for other in later:
        days = other.mid - seed.mid
        if not .7 <= days <= 6.5 or other.night == seed.night:
            continue
        # Fast prefilter on the observed sky directions; injected objects pass the identical gate.
        if np.linalg.norm(other.sky - unit(seed.sky + days*seed.rate)) > np.deg2rad(.25 + .15*days**2):
            continue
        pred = unit(short_arc(pos,vel,days)-other.site)
        error = np.linalg.norm(pred - other.sky,axis=1)/ARC
        if np.min(error) < 600:
            links.append((other, error))
    for (b, e1), (c, e2) in combinations(links, 2):
        if not .7 <= c.mid-b.mid <= 3.25 or b.mid-seed.mid > 3.25 or b.night == c.night:
            continue
        i = np.argmin(e1**2+e2**2)
        if max(e1[i],e2[i]) < 600:
            yield (seed,b,c), np.r_[pos[i], vel[i]]


def separate(x, y):
    return np.linalg.norm(x-y,axis=1)/ARC


def correct(tracks, guess):
    seed = tracks[0]
    times = np.concatenate([a.t for a in tracks]) - seed.mid
    obs = np.concatenate([a.u for a in tracks])
    sites = np.array([observer(t+seed.mid) for t in times])
    def calc(state):
        # State epoch is the seed reception time minus its own range/c.
        delay0 = np.linalg.norm(state[:3]-seed.site)/C
        pos = short_arc(state[:3],state[3:],times[:,None])
        dates = times + delay0 - np.linalg.norm(pos-sites,axis=1)/C
        return unit(motion(state,dates)-sites)
    def resid(state):
        return ((calc(state)-obs)/ARC).ravel()
    try:
        fit = least_squares(resid,guess,x_scale=np.array([1.]*3+[.03]*3),
                            ftol=1e-7,xtol=1e-9,gtol=1e-6,max_nfev=70)
        sep = separate(calc(fit.x),obs)
        r,v = fit.x[:3], fit.x[3:]
        excess = (v@v - 2*MU/np.linalg.norm(r)) / KMS**2
        rms = float(np.sqrt(np.mean(sep**2)))
        if fit.success and rms <= 1.0 and max(sep) <= 2.5 and excess > 100:
            return dict(ids=[a.name for a in tracks], rms_arcsec=rms,
                        max_arcsec=float(max(sep)), v_inf_helio_kms=float(np.sqrt(excess)),
                        distance_au=float(np.linalg.norm(r-seed.site)),
                        seed_utc_days_past_j2000=float(seed.mid), state_au_au_day=fit.x.tolist())
    except (ValueError, FloatingPointError, OverflowError):
        return None


def scan(tracks):
    result, n = [], 0
    tracks = sorted(tracks, key=lambda a:a.mid)
    tt = np.array([a.mid for a in tracks])
    for i,a in enumerate(tracks):
        for triple,guess in proposals(a, tracks[i+1:np.searchsorted(tt,a.mid+6.5,side="right")]):
            n += 1
            found = correct(triple,guess)
            if found is not None:
                result.append(found)
    return n,result



def make_injection(rng, tracks, idx, epoch_triples):
    # Force three tracklets at real timestamps. Sky at the first epoch is real,
    # later positions need not be inside an acquired exposure.
    a,b,c = epoch_triples[rng.integers(len(epoch_triples))]
    obs0 = a.site
    rho = np.exp(rng.uniform(np.log(.45),np.log(2.8)))
    position = obs0 + rho*a.sky
    vinf = rng.uniform(15.,60.)
    direction = unit(rng.normal(size=3))
    velocity = direction*np.sqrt((vinf*KMS)**2 + 2*MU/np.linalg.norm(position))
    state = np.r_[position, velocity]
    parts = []
    for j,tr in enumerate([a,b,c]):
        off = tr.t - a.mid
        sites = np.array([observer(float(t)) for t in tr.t])
        dates = off.copy()
        for _ in range(3):
            pos = motion(state,dates)
            dates = off+(rho-np.linalg.norm(pos-sites,axis=1))/C
        pred = unit(motion(state,dates)-sites)
        noisy = unit(pred+rng.normal(0,.2*ARC,pred.shape))
        if SURVEY=="F51":
            ra=np.rad2deg(np.arctan2(noisy[:,1],noisy[:,0]))%360
            dec=np.rad2deg(np.arcsin(np.clip(noisy[:,2],-1,1)))
            keep=(40<=ra)&(ra<50)&(10<=dec)&(dec<20)
            if np.sum(keep)<2:
                parts.append(None)
                continue
            parts.append(attributable(f"inj-{idx}-{j}",tr.t[keep],noisy[keep]))
        else:
            parts.append(attributable(f"inj-{idx}-{j}",tr.t,noisy))
    return parts,dict(rho_au=float(rho), true_vinf_kms=float(vinf), times_utc_days=[float(t.mid) for t in [a,b,c]])


def cadence_templates(tracks):
    # Sample each eligible triple of *nights* equally, then a real track on each.
    nights = collections.defaultdict(list)
    for a in tracks:
        nights[a.night].append(a)
    valid=[]
    for a,b,c in combinations(sorted(nights),3):
        if b-a > 4 or c-b > 4:
            continue
        valid.append((a,b,c))
    return nights,valid


def draw_times(rng, by_day, nights):
    options=[by_day[i] for i in nights[rng.integers(len(nights))]]
    return [x[rng.integers(len(x))] for x in options]


def rotate_sky(tracks,rng):
    angles = {a.night:rng.uniform(-3,3)*np.pi/180 for a in sorted(tracks,key=lambda a:a.night)}
    copies=[]
    for a in tracks:
        angle=angles[a.night]
        co,si=np.cos(angle),np.sin(angle)
        rot=np.array([[co,-si,0],[si,co,0],[0,0,1]])
        copies.append(attributable(a.name, a.t, a.u@rot.T))
    return copies


def main():
    global SURVEY
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survey",choices=sorted(SITES),default="X05")
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--trials",type=int,default=100)
    args=parser.parse_args()
    SURVEY=args.survey
    path=args.data/f"{SURVEY.lower()}.mpc"
    real,qc=parse(path)
    before,accepted=scan(real)
    rng=np.random.default_rng(RNG_SEED)
    by_day,options=cadence_templates(real)
    # Remove infeasible time gaps at the epoch level before drawing an object.
    templates=[]
    if options:
        for _ in range(1000):
            trio=draw_times(rng,by_day,options)
            if all(.7<=y.mid-x.mid<=3.25 for x,y in zip(trio,trio[1:])):
                templates.append(trio)
    result=dict(survey=SURVEY,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),qc=qc,
                real_eligible_per_utc_night={str(J2000+dt.timedelta(days=int(k))):len(v) for k,v in sorted(by_day.items())},
                real_3night_proposals=before,real_passing=accepted,seed=RNG_SEED)
    recovery=[]
    for idx in range(args.trials if templates else 0):
        injected,row=make_injection(rng,real,idx,templates)
        row["eligible"]=all(a is not None for a in injected)
        row["recovered"]=False
        row["mixed_accepts"]=0
        if row["eligible"]:
            # The linker considers every qualifying triple independently (no
            # competition). Restrict repeated work to triples touching an injection.
            pool=sorted(real+injected,key=lambda a:a.mid)
            inject_ids={x.name for x in injected}
            for i,a in enumerate(pool):
                if a.name not in inject_ids and not any(.7<=z.mid-a.mid<=6.5 for z in injected):
                    continue
                for triple,guess in proposals(a,pool[i+1:]):
                    ids={x.name for x in triple}
                    if ids.isdisjoint(inject_ids):
                        continue
                    found=correct(triple,guess)
                    if found and ids == inject_ids:
                        row["recovered"]=True
                    elif found:
                        row["mixed_accepts"]+=1
        recovery.append(row)
    n=len(recovery)
    k=sum(x["recovered"] for x in recovery)
    result.update(injected=n,eligible_injections=sum(x["eligible"] for x in recovery),
                  recovered=k,mixed_passes=sum(x["mixed_accepts"] for x in recovery),
                  exact_95pct_recovery_interval=([0. if k==0 else float(beta.ppf(.025,k,n-k+1)),
                                                  1. if k==n else float(beta.ppf(.975,k+1,n-k))] if n else None))
    scramble=[]
    for _ in range(5):
        nprop,hits=scan(rotate_sky(real,rng))
        scramble.append([nprop,len(hits)])
    result["night_rotation_controls_proposed_passing"]=scramble
    (args.data/f"{SURVEY.lower()}-injections.json").write_text(json.dumps(recovery,indent=2)+"\n")
    (args.data/f"{SURVEY.lower()}-result.json").write_text(json.dumps(result,indent=2)+"\n")
    print("RESULT_JSON="+json.dumps(result),flush=True)

if __name__=="__main__":
    main()
