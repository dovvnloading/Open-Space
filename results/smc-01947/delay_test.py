"""Apply the published eclipse-delay phase to an independent brightness clock."""
import json

import numpy as np
from scipy.optimize import least_squares

from measure import EPOCH, PERIOD, orbit, series

FREQUENCY = 5.644829612245024
HONG_EPOCH = 5002.2972
HONG_BINARY_PERIOD = 2.401342


def estimate(t, magnitude, error, rate=.34, phase_degrees=192):
    # Omit broad windows around both eclipses; apsidal motion changes their shapes.
    phase = (t-10001.8882)/PERIOD % 1
    use = (abs((phase+.5) % 1-.5) > .085) & (abs(phase-.5) > .15)
    t, magnitude, error = t[use], magnitude[use], error[use]
    nuisance = orbit(t,3)
    cycle = (t-HONG_EPOCH)/HONG_BINARY_PERIOD
    predicted = np.sin(np.deg2rad(rate*cycle+phase_degrees))
    sigma = np.hypot(error,.004)

    def model_and_jac(c):
        freq = FREQUENCY+c[0]*1e-5
        elapsed = t-EPOCH-c[1]*predicted/1440
        angle = 2*np.pi*freq*elapsed
        co, si = np.cos(angle), np.sin(angle)
        dangle = -c[-2]*si+c[-1]*co
        y = c[2]+nuisance@c[3:-2]+c[-2]*co+c[-1]*si
        jac = np.column_stack((dangle*2*np.pi*elapsed*1e-5,
                               -dangle*2*np.pi*freq*predicted/1440,
                               np.ones(len(t)),nuisance,co,si))
        return y, jac

    start = np.r_[0,0,0,np.zeros(6),.008,-.003]
    fit = least_squares(lambda c:(model_and_jac(c)[0]-magnitude)/sigma,start,
                        jac=lambda c:model_and_jac(c)[1]/sigma[:,None],
                        loss='soft_l1',f_scale=2)
    if not fit.success:
        raise RuntimeError(fit.message)
    val, jac = model_and_jac(fit.x)
    jac /= sigma[:,None]
    r = (magnitude-val)/sigma
    weight = 1+(r/2)**2
    influence = jac*(r/np.sqrt(weight))[:,None]
    bread = np.linalg.pinv(jac.T@(jac*weight[:,None]**-1.5))
    nights,group = np.unique(np.floor(t),return_inverse=True)
    score = np.zeros((len(nights),len(start)))
    np.add.at(score,group,influence)
    se = np.sqrt(np.sum((score @ bread[:,1])**2))
    return {'observations':len(t),'nights':len(nights),'delay_minutes':float(fit.x[1]),
            'night_cluster_se_minutes':float(se)}


def analyze():
    results = {}
    grid = []
    for label,t,y,e in series():
        if label not in ('OGLE-4-I','OGLE-3-I','EROS-sm0010n21934-R','EROS-sm0010n21934-B'):
            continue
        results[label] = estimate(t,y,e)
        if label == 'OGLE-4-I':
            for rate in [.32,.34,.36]:
                for phase in [163,192,221]:
                    result = estimate(t,y,e,rate,phase)
                    grid.append({'rate_deg_per_cycle':rate,'phase_deg':phase,**result})
    return {'published_central_delay_minutes':.008*1440,'fits':results,'ogle4_sensitivity':grid}


if __name__ == '__main__':
    print(json.dumps(analyze(),indent=2))
