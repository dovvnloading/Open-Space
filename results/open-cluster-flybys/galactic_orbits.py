"""Check the nominal cluster orbits in galpy MWPotential2014 (astropy and galpy)."""
import json
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

from replay import HERE, INPUTS, normalise, screen


def trajectory(phase, times):
    a,b,d,ma,md,rv = phase
    orbit = Orbit([a*u.deg,b*u.deg,d*u.pc,ma*u.mas/u.yr,md*u.mas/u.yr,rv*u.km/u.s],
                  radec=True,ro=8,vo=220,zo=.0208,solarmotion='schoenrich')
    orbit.integrate(times*u.Myr, MWPotential2014, method='dop853')
    return np.array([coordinate(times*u.Myr) for coordinate in (orbit.x,orbit.y,orbit.z)]).T * 1000


def calculate():
    hunt=pd.read_csv(INPUTS/'hunt_reffert.csv')
    ucc=pd.read_csv(INPUTS/'ucc_centres.csv').set_index('name')
    candidates=screen(hunt)
    hunt=hunt.set_index('Name')
    checks=[]
    for row in candidates.itertuples(index=False):
        times=np.linspace(0,8*np.sign(row.t_myr),401)
        for label in ['Hunt-Reffert','UCC']:
            phases=[]
            for name in [row.first,row.second]:
                if label=='UCC':
                    g=ucc.loc[normalise(name)]
                    phases.append([g.RA_ICRS,g.DE_ICRS,1000/(g.Plx+.017),g.pmRA,g.pmDE,g.Rv])
                else:
                    g=hunt.loc[name]
                    phases.append([g.RAdeg,g.DEdeg,g.dist50,g.pmRA,g.pmDE,g.RV])
            delta=trajectory(phases[0],times)-trajectory(phases[1],times)
            order=np.argsort(times)
            curve=CubicSpline(times[order],delta[order])
            distance=lambda t:float(np.linalg.norm(curve(t)))
            opt=minimize_scalar(distance,bounds=(min(times),max(times)),method='bounded')
            checks.append(dict(first=row.first,second=row.second,catalog=label,
                               pericentre_pc=round(float(opt.fun),4),t_myr=round(float(opt.x),4),
                               grid_min_pc=round(min(distance(t) for t in times[::100]),4)))
    return checks


if __name__=='__main__':
    print(json.dumps(calculate(),indent=2))
