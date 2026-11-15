"""Replicate the source and receiver epochs with three dynamical calculations."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, root

from orbits import AU_KM, INPUTS, K, catalog, date_string, day, frozen_jpl, state_spline


def miriade(name):
    rows = json.loads((INPUTS/f"miriade_{name}.json").read_text())["data"]
    keys = ["Date","px","py","pz","vx","vy","vz"]
    return state_spline(np.array([[r[k] for k in keys] for r in rows]))


def mpc_states():
    data = np.loadtxt(INPUTS/"planets.csv", delimiter=",", skiprows=1)
    tracks = np.array([data[data[:,0]==i,1:] for i in range(1,9)])
    planets = CubicHermiteSpline(tracks[0,:,0]-2451545, tracks[:,:,1:4].swapaxes(0,1),
                                 tracks[:,:,4:7].swapaxes(0,1),extrapolate=False)
    mass_ratios = np.array([6023600,408523.71,328900.56,3098708,1047.3486,3497.898,22902.98,19412.24])
    mu = K**2/mass_ratios
    def rhs(t, x):
        r, v = x[:3], x[3:]
        bodies = planets(t)
        dr = bodies-r
        acc = -K**2*r/np.linalg.norm(r)**3
        acc += np.sum(mu[:,None]*(dr/np.linalg.norm(dr,axis=1)[:,None]**3 -
                                  bodies/np.linalg.norm(bodies,axis=1)[:,None]**3),axis=0)
        return np.r_[v,acc]

    names,tp,q,e,epochs,P,Q = catalog()
    splines = []
    for name in ("3I/ATLAS","383P/Christensen"):
        i = np.flatnonzero(names==name).item()
        ecc, epoch = e[i], epochs[i]
        a = q[i]/abs(1-ecc)
        mean = K*(epoch-tp[i])/a**1.5
        if ecc>1:
            h = brentq(lambda z:ecc*np.sinh(z)-z-mean,-10,10)
            b = np.sqrt(ecc**2-1)
            x,y = a*(ecc-np.cosh(h)),a*b*np.sinh(h)
            factor = K/np.sqrt(a)/(ecc*np.cosh(h)-1)
            vx,vy = -factor*np.sinh(h),factor*b*np.cosh(h)
        else:
            m = (mean+np.pi)%(2*np.pi)-np.pi
            h = brentq(lambda z:z-ecc*np.sin(z)-m,-np.pi,np.pi)
            b = np.sqrt(1-ecc**2)
            x,y = a*(np.cos(h)-ecc),a*b*np.sin(h)
            factor = K/np.sqrt(a)/(1-ecc*np.cos(h))
            vx,vy = -factor*np.sin(h),factor*b*np.cos(h)
        y0 = np.r_[x*P[i]+y*Q[i],vx*P[i]+vy*Q[i]]
        sol = solve_ivp(rhs,(epoch,day("2025-10-01")),y0,rtol=1e-12,atol=1e-14,
                        method="DOP853",dense_output=True,max_step=4)
        if not sol.success:
            raise RuntimeError(sol.message)
        splines.append(lambda t,s=sol.sol:s(t)[:3])
    return splines


def crossing(source, receiver):
    seed = source(day("2025-10-04"))
    u = seed/np.linalg.norm(seed)
    east = np.cross([0,0,1],u)
    east /= np.linalg.norm(east)
    axes = np.array([east,np.cross(u,east)])
    def difference(times):
        a,b = source(times[0]),receiver(times[1])
        return axes@(a/np.linalg.norm(a)-b/np.linalg.norm(b))
    fit = root(difference,[day("2025-10-04"),day("2025-10-11")],tol=1e-12)
    t1,t2 = fit.x
    a,b = source(t1),receiver(t2)
    if np.linalg.norm(difference(fit.x))>1e-10 or a@b<=0 or t1>=t2:
        raise RuntimeError("Directional root failed")
    distance = np.linalg.norm(b)-np.linalg.norm(a)
    speed = distance*AU_KM/(86400*(t2-t1))
    return t1,t2,distance,speed


if __name__ == "__main__":
    for name, pair in [("Horizons",[frozen_jpl(x) for x in ("3I","383P")]),
                       ("MPC gravity",mpc_states()),("Miriade",[miriade(x) for x in ("3I","383P")])]:
        t1,t2,distance,speed = crossing(*pair)
        print(name,date_string(t1),date_string(t2),f"{distance:.8f} au {speed:.5f} km/s")
        if name=="Horizons":
            a,b = pair[0](t1),frozen_jpl("499")(t1)
            sep = np.rad2deg(np.arccos(a@b/(np.linalg.norm(a)*np.linalg.norm(b))))
            print(f"Mars proxy is {sep:.4f} deg, {np.linalg.norm(a-b):.6f} au from 3I")
