"""Select exactly those circular grid rows that can match an injected orbit."""
import numpy as np
from astropy.time import Time
from scipy.spatial import cKDTree
from geometry import AU,C,ARCSEC,basis,frequency,tangent_basis


def compatible_rows(cube,plants,epoch,step=500.,phase_step=1.,match=2.):
    # The first filter alone gives a superset. The full test later enforces both.
    j=cube['path_mjd'].shape[1]//2
    ra,dec,delta=(cube[k][:,j] for k in ['path_ra','path_dec','path_delta'])
    times=Time(cube['path_mjd'][:,j],format='mjd',scale='utc').tdb.jd-delta*AU/C/86400
    los,east,north=tangent_basis(ra,dec);u,v,_=basis(str(cube['planet']),epoch)
    alpha=np.arange(0.,360.,phase_step)
    tree=cKDTree(np.column_stack((plants['x'][:,0],plants['y'][:,0])))
    rows=[]
    for radius in np.arange(30000.,130000.+1.,step):
        an=np.deg2rad(alpha[:,None])+frequency(str(cube['planet']),radius)*(times-epoch)*86400
        pos=radius*(np.cos(an)[:,:,None]*u+np.sin(an)[:,:,None]*v)
        denom=delta*AU+(pos*los).sum(-1)
        x=cube['center_xy'][0]-(pos*east).sum(-1)/denom*ARCSEC/float(cube['scale_arcsec'])
        y=cube['center_xy'][1]+(pos*north).sum(-1)/denom*ARCSEC/float(cube['scale_arcsec'])
        for k,possible in enumerate(tree.query_ball_point(np.column_stack((x[:,0],y[:,0])),match+1e-3)):
            if possible and np.any(np.hypot(x[k]-plants['x'][possible],y[k]-plants['y'][possible]).max(1)<=match+.001):
                rows.append([radius,0.,0.,0.,0.,alpha[k]])
    return rows
