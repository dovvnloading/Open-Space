"""Compare 60/100-micron SCANPI timelines at stationary and moving positions."""
import argparse
import csv
from datetime import datetime
import importlib.util
import json
from pathlib import Path
import re
import tarfile
import io

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np

ROOT = Path(__file__).resolve().parent
CENTER = SkyCoord(35.74075,-48.5125,unit='deg')
GALAXY = SkyCoord(35.7544007,-48.5112676,unit='deg')


def offsets(xy, pa):
    dx,dy = CENTER.spherical_offsets_to(xy)
    th=np.deg2rad(pa)
    return dx.arcmin*np.sin(th)+dy.arcmin*np.cos(th),dx.arcmin*np.cos(th)-dy.arcmin*np.sin(th)


def filtered(x,f,width,window_fwhm=3):
    """Fixed-width template amplitude with local constant and slope nuisance terms.

    The symmetric window makes the linear term orthogonal to the template.
    Correlated oversampling is carried into off-source measurements.
    """
    step=np.median(np.diff(x))
    assert np.allclose(np.diff(x),step)
    n=int(round(window_fwhm*width/step))
    dx=np.arange(-n,n+1)*step
    kernel=np.exp(-4*np.log(2)*(dx/width)**2)
    kernel-=kernel.mean()
    return np.convolve(f,kernel/(kernel@kernel),'same')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--orbit-dir',required=True,type=Path)
    args=ap.parse_args()
    s=importlib.util.spec_from_file_location('orbit',args.orbit_dir/'propagate.py')
    orbit=importlib.util.module_from_spec(s)
    s.loader.exec_module(orbit)
    states=orbit.family(json.loads((args.orbit_dir/'orbit_input.json').read_text()))
    july=orbit.sky(states,orbit.epoch('1983-07-09T23:50:30'))
    anchor=offsets(SkyCoord(*july,unit='deg'),332.)[0]
    gate=abs(anchor)<.5  # descriptive July-compatible slice
    groupout=[]
    summary=json.loads((ROOT/'scanpi'/'sign.json').read_text()) if (ROOT/'scanpi'/'sign.json').exists() else {}
    with tarfile.open(ROOT/'scanpi'/'profiles.tar.gz') as archive:
        for band,width in [(3,1.5),(4,3.)]:
            tab=Table.read(ROOT/'scanpi'/f'b{band}_summary.tbl',format='ascii.ipac')
            results=[]
            for r in tab:
                if r['scan']>=999:continue
                n=int(r['scan'])
                member=f'b{band}/'+Path(re.search(r'href="([^"]+)"',r['table']).group(1)).name
                raw=archive.extractfile(member).read().decode()
                data=np.array([[float(z) for z in line.split()[:2]] for line in raw.splitlines() if line.startswith(' ')])
                x,signal=data.T
                amp=filtered(x,signal,width)
                pa=float(re.search(r'^\\pa\s*=\s*([-\d.]+)',raw,re.M).group(1))
                date=datetime.strptime(r['utcs_start'],'%m/%d/%y  %H:%M:%S.%f')
                ra,dec=orbit.sky(states,orbit.epoch(date.isoformat()))
                along,cross=offsets(SkyCoord(ra,dec,unit='deg'),pa)
                gal,galcross=offsets(GALAXY,pa)
                controls=np.r_[np.arange(-40.,-8.,4.),np.arange(12.,43.,4.)]
                field=np.interp(controls,x,amp)
                sig=1.4826*np.median(abs(field-np.median(field)))
                out=dict(band_um={3:60,4:100}[band],scan=n,utc=date.isoformat(),pa_deg=pa,
                         scan_offset_arcmin=r['offset'],archive_peak_noise=r['snr'],
                         archive_miss_arcmin=r['miss'],archive_template_Jy=r['amp'],
                         fixed_Jy=np.interp(0,x,amp),control_median_Jy=np.median(field),
                         control_robust_scatter_Jy=sig,gal_along_arcmin=gal,gal_cross_arcmin=galcross,
                         july_gated_along_min_arcmin=along[gate].min(),july_gated_along_max_arcmin=along[gate].max(),
                         july_gated_cross_min_arcmin=cross[gate].min(),july_gated_cross_max_arcmin=cross[gate].max())
                for sign in [-1,1]:
                    y=np.interp(sign*along[gate],x,amp)
                    out.update({f'gal_Jy_sign{sign:+}':np.interp(sign*gal,x,amp),
                                f'moving_Jy_min_sign{sign:+}':y.min(),f'moving_Jy_median_sign{sign:+}':np.median(y),
                                f'moving_Jy_max_sign{sign:+}':y.max()})
                group=('June' if date.month==6 else 'July' if date.month==7 else 'Aug_Sep')
                results.append((out,group,x,amp,along[gate],gal,sig))
                groupout.append(out)
            for group in ['June','July','Aug_Sep']:
                same=[r for r in results if r[1]==group]
                controls=np.r_[np.arange(-40.,-8.,4.),np.arange(12.,43.,4.)]
                def avg(xx):
                    return np.average(xx,axis=0,weights=[1/max(r[6],.04)**2 for r in same])
                const=avg([np.interp(0,r[2],r[3]) for r in same])
                blank=avg([np.interp(controls,r[2],r[3]) for r in same])
                for sign in [-1,1]:
                    moving=avg([np.interp(sign*r[4],r[2],r[3]) for r in same])
                    g=avg([np.interp(sign*r[5],r[2],r[3]) for r in same])
                    summary[f'{band}_{group}_sign{sign:+}']=dict(n=len(same),fixed_Jy=const,galaxy_Jy=g,
                        moving_min_med_max_Jy=np.percentile(moving,[0,50,100]).tolist(),
                        blank_median_Jy=np.median(blank),blank_MAD_sigma_Jy=1.4826*np.median(abs(blank-np.median(blank))),
                        gated_models=int(gate.sum()),model_count=len(gate))
    with (ROOT/'scanpi'/'track_measurements.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,groupout[0]);w.writeheader();w.writerows(groupout)
    (ROOT/'scanpi'/'track_summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
