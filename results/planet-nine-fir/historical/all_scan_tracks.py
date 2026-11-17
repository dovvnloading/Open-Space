"""Check all nominal motion tracks against both independent temporal groups at 60 um."""
import argparse
from datetime import datetime
import importlib.util
import json
from pathlib import Path
import re
import tarfile

from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
from scan_tracks import ROOT, offsets, filtered


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orbit-dir',type=Path,required=True)
    parser.add_argument('--width',type=float,default=1.5)
    parser.add_argument('--window',type=float,default=3.)
    parser.add_argument('--output',type=Path,default=ROOT/'scanpi/all_scan_tracks.json')
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('orbit',args.orbit_dir/'propagate.py')
    orbit=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orbit)
    states=orbit.family(json.loads((args.orbit_dir/'orbit_input.json').read_text()))
    tab=Table.read(ROOT/'scanpi/b3_summary.tbl',format='ascii.ipac')
    displacements=np.r_[np.arange(-40.,-8.,4.),np.arange(12.,43.,4.)]
    groups={'july':[],'late':[]}
    with tarfile.open(ROOT/'scanpi/profiles.tar.gz') as archive:
        for r in tab:
            if r['scan']>=999:continue
            date=datetime.strptime(r['utcs_start'],'%m/%d/%y  %H:%M:%S.%f')
            if date.month<7:continue
            group='july' if date.month==7 else 'late'
            member='b3/'+Path(re.search(r'href="([^"]+)"',r['table']).group(1)).name
            raw=archive.extractfile(member).read().decode()
            data=np.array([[float(z) for z in line.split()[:2]] for line in raw.splitlines() if line.startswith(' ')])
            x,signal=data.T
            amp=filtered(x,signal,args.width,args.window)
            pa=float(re.search(r'^\\pa\s*=\s*([-\d.]+)',raw,re.M).group(1))
            ra,dec=orbit.sky(states,orbit.epoch(date.isoformat()))
            along,cross=offsets(SkyCoord(ra,dec,unit='deg'),pa)
            blank=np.interp(displacements,x,amp)
            sigma=1.4826*np.median(abs(blank-np.median(blank)))
            row=np.interp(along,x,amp)
            groups[group].append((1/max(sigma,.04)**2, np.interp(0,x,amp), row, blank))
    out={}
    for group,samples in groups.items():
        w,f,m,b=zip(*samples)
        flux=np.average(m,axis=0,weights=w)
        controls=np.average(b,axis=0,weights=w)
        center=np.median(controls)
        scatter=1.4826*np.median(abs(controls-center))
        out[group]=(flux,center,scatter,np.average(f,weights=w))
    # A persistent equal-brightness planet should give positive amplitudes in both intervals.
    july,late=(out[g] for g in ('july','late'))
    minimum=np.minimum(july[0],late[0])
    iz=int(np.argmax(minimum))
    def item(i):return {'grid_index':int(i),'July_Jy':july[0][i],'late_Jy':late[0][i]}
    result={'fwhm_arcmin':args.width,'window_halfwidth_arcmin':args.width*args.window,'grid_size':len(minimum),'largest_minimum_flux_Jy':float(max(minimum)),
            'best_balanced_unweighted':item(iz),
            'grid_meeting_flux_threshold_in_both_groups':{
              str(t):int(np.sum(minimum>=t)) for t in (.05,.10,.12,.15,.18,.20,.22,.24)},
            'grid_July_flux_at_least_0p20_Jy':int(np.sum(july[0]>=.20)),
            'max_late_Jy_at_July_at_least_0p20':max(late[0][july[0]>=.20]),
            'fixed_July_late_Jy':[float(july[3]),float(late[3])]}
    args.output.write_text(json.dumps(result,indent=2)+'\n')



if __name__=='__main__':main()
