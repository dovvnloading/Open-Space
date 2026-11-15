"""Re-run the close-approach screen for both recent candidate eclipse dates."""
from io import StringIO
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar

DAY = 86400.0


def interpolator(body, year, step):
    params=dict(COMMAND=body,CENTER='500@899',EPHEM_TYPE='VECTORS',
                START_TIME=f'{year}-01-01',STOP_TIME=f'{year+20}-01-05',
                STEP_SIZE=step,OUT_UNITS='KM-S',REF_PLANE='FRAME',TIME_TYPE='TDB',
                VEC_TABLE='2',VEC_CORR='NONE',CSV_FORMAT='YES',OBJ_DATA='NO')
    args={k:"'"+v+"'" for k,v in params.items()}
    url='https://ssd.jpl.nasa.gov/api/horizons.api?format=text&'+urlencode(args)
    text=urlopen(url,timeout=120).read().decode()
    rows=text.split('$$SOE')[1].split('$$EOE')[0]
    a=np.loadtxt(StringIO(rows),delimiter=',',usecols=(0,2,3,4,5,6,7))
    return CubicHermiteSpline(a[:,0],a[:,1:4],a[:,4:]*DAY)


def run(year):
    triton=interpolator('801',year,'1d')
    nereid=interpolator('802',year,'4d')
    grid=np.arange(triton.x[0],triton.x[-1]-4,1/16)
    for source in ('10','399'):
        observer=interpolator(source,year,'4d')
        def miss(t):
            moon=nereid(t)
            ray=observer(t)-moon
            n=ray/np.linalg.norm(ray,axis=-1,keepdims=True)
            return np.linalg.norm(np.cross(triton(t)-moon,n),axis=-1)
        d=miss(grid)
        indices=np.flatnonzero((d[1:-1]<d[:-2])&(d[1:-1]<d[2:])&(d[1:-1]<35000))+1
        for i in indices:
            opt=minimize_scalar(lambda h:miss(grid[i]+h),bounds=(-.1,.1),method='bounded')
            if opt.fun<12000:
                jd=grid[i]+opt.x
                t=np.datetime64('2000-01-01T12:00:00')+np.timedelta64(round((jd-2451545)*DAY),'s')
                print(source,t,round(float(opt.fun)))

if __name__=='__main__':
    run(2000)
    run(2020)
