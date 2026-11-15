"""Point-mass examples for HD 73712 in au, solar mass units and years."""
from itertools import product
from concurrent.futures import ProcessPoolExecutor
import json
import math
import numpy as np
from scipy.integrate import solve_ivp

G = 4 * math.pi ** 2
M = np.array([1 / 1.86, .86 / 1.86, .7])
PI = 48.717 / 365.25
def outer_period(e):
    return 20.96 if e == .98 else 20.83

def initial(eout, inc, phase):
    # Both osculating Jacobi ellipses are specified at outer apastron.
    ai = (PI ** 2 * sum(M[:2])) ** (1 / 3)
    ao = (outer_period(eout) ** 2 * sum(M)) ** (1 / 3)
    ei = .11
    f = math.radians(phase)
    co, si = math.cos(f), math.sin(f)
    ci, sii = math.cos(math.radians(inc)), math.sin(math.radians(inc))
    inner_r = ai * (1 - ei ** 2) / (1 + ei * co) * np.array([co, si*ci, si*sii])
    inner_v = math.sqrt(G * sum(M[:2]) / (ai*(1-ei**2))) * np.array([-si, (co+ei)*ci, (co+ei)*sii])
    outer_r = np.array([-ao*(1+eout), 0., 0.])
    outer_v = np.array([0., -math.sqrt(G*sum(M)/ao * (1-eout)/(1+eout)), 0.])
    def to_bary(inner, outer):
        inner_cm = -M[2] / sum(M) * outer
        return np.array([inner_cm-M[1]*inner, inner_cm+M[0]*inner, inner_cm+outer])
    return np.concatenate([to_bary(inner_r, outer_r), to_bary(inner_v, outer_v)]).ravel()

def acceleration(t, state):
    x0,y0,z0,x1,y1,z1,x2,y2,z2, *vel = state
    dx1,dy1,dz1 = x1-x0,y1-y0,z1-z0
    dx2,dy2,dz2 = x2-x0,y2-y0,z2-z0
    dx3,dy3,dz3 = x2-x1,y2-y1,z2-z1
    c1 = G/(dx1*dx1+dy1*dy1+dz1*dz1)**1.5
    c2 = G/(dx2*dx2+dy2*dy2+dz2*dz2)**1.5
    c3 = G/(dx3*dx3+dy3*dy3+dz3*dz3)**1.5
    fac1,fac2 = c1*M[1],c2*M[2]
    fac3,fac4 = -c1*M[0],c3*M[2]
    fac5,fac6 = -c2*M[0],-c3*M[1]
    return [*vel,fac1*dx1+fac2*dx2,fac1*dy1+fac2*dy2,fac1*dz1+fac2*dz2,
            fac3*dx1+fac4*dx3,fac3*dy1+fac4*dy3,fac3*dz1+fac4*dz3,
            fac5*dx2+fac6*dx3,fac5*dy2+fac6*dy3,fac5*dz2+fac6*dz3]

def energy(y):
    r,v = y[:9].reshape(3,3),y[9:].reshape(3,3)
    return .5*np.sum(M[:,None]*v*v) - sum(G*M[i]*M[j]/np.linalg.norm(r[i]-r[j]) for i,j in [(0,1),(0,2),(1,2)])

def inner_period(y):
    r, v = y[:9].reshape(3, 3), y[9:].reshape(3, 3)
    delta = .5*np.sum((v[1]-v[0])**2) - G/np.linalg.norm(r[1]-r[0])
    if delta >= 0:
        return None
    return 365.25*(-G/(2*delta))**1.5


def run(args):
    e, inc, phase, rtol = args
    y = initial(e, inc, phase)
    # Evaluate after first periastron. This is an osculating two-body period.
    t_end = .72 * outer_period(e)
    sol = solve_ivp(acceleration, (0., t_end), y, method="DOP853",
                    rtol=rtol, atol=rtol/100, t_eval=[t_end])
    row = dict(e=e, mutual_inclination_deg=inc, true_anomaly_deg=phase, rtol=rtol,
               integration_success=sol.success)
    if sol.success:
        final = sol.y[:, -1]
        row.update(inner_period_days=inner_period(final), fractional_energy_error=float(abs(energy(final)/energy(y)-1)))
    return row


if __name__ == "__main__":
    scenarios = product([.98,.93,.75], [0,90,180], [0,90,180,270], [1e-10,1e-12])
    with ProcessPoolExecutor(4) as pool:
        for row in pool.map(run, scenarios):
            print(json.dumps(row), flush=True)
