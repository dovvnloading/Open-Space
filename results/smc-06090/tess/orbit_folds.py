"""Phase-bin post-detrending aperture data within each orbital cycle for timing."""
from pathlib import Path
import inspect
import json
import numpy as np
import photometry as phot

# Use the same callable as the pixel stack. Its optional return retains native
# photometry in memory; the exported files contain only sums and counts per orbit.
def run():
    store={}
    for p in sorted((phot.HERE/'raw').glob('*.fits')):
        t,ph,aps,orb=phot.lightcurves(p,return_native=True)
        night=np.unique(orb)
        S=np.zeros((len(night),120,3));N=np.zeros((len(night),120))
        b=np.floor((ph%1)*120).astype(int)
        for i,cycle in enumerate(night):
            sel=orb==cycle
            N[i]=np.bincount(b[sel],minlength=120)
            for k in range(3):S[i,:,k]=np.bincount(b[sel],weights=aps[sel,k],minlength=120)
        f=p.stem.lstrip('s0')
        for key,val in [('sums',S),('counts',N),('cycles',night)]:store['s'+f+'_'+key]=val
    np.savez_compressed(phot.HERE/'orbital_folds.npz',**store)

if __name__=='__main__':run()
