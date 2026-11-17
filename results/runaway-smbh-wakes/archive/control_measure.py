"""Report the targeted positive control and ground-based forced-aperture diagnostic."""
import argparse,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from wake_quality import extended_profile

HOST=SkyCoord('02h41m45.43s -08d20m55.4s')
PA=147*u.deg  # SSE side of the published 327-degree slit axis
TIP=HOST.directional_offset_by(PA,7.8*u.arcsec)

def measure(images,results):
    report={'host_deg':[HOST.ra.deg,HOST.dec.deg],'pa_east_of_north_deg':147,'ACS':[], 'Legacy':[]}
    for band,obs in [('606','jety02010'),('814','jety02020')]:
        rows=json.loads((results/f'control_{band}.json').read_text())['detections']
        with fits.open(images/f'{obs}_drc.fits') as hdul:
            w=WCS(hdul['SCI'].header)
            sl=np.s_[1320:1704,2200:2584]
            science=hdul['SCI'].data[sl].copy(); valid=hdul['WHT'].data[sl]>0
        host_xy=np.array(w.world_to_pixel(HOST));tip_xy=np.array(w.world_to_pixel(TIP))
        tangent=tip_xy-host_xy;tangent/=np.linalg.norm(tangent)
        image_pa=np.degrees(np.arctan2(tangent[1],tangent[0]))%180
        for rank,row in enumerate(rows,1):
            pos=np.array([row['x'],row['y']])-host_xy
            s=pos@tangent*.05
            cross=abs(pos@[-tangent[1],tangent[0]]*.05)
            angle=abs((row['pa_image_deg']-image_pa+90)%180-90)
            if 1 <= s <= 7 and cross < .3 and angle <= 8:
                local=dict(row)
                for k in ('x0','x1'):local[k]-=2200
                for k in ('y0','y1'):local[k]-=1320
                dist=min(np.hypot(row['x'+str(i)]-host_xy[0],row['y'+str(i)]-host_xy[1]) for i in (0,1))*.05
                report['ACS'].append(dict(band=band,rank=rank,total=len(rows),axial_distance_arcsec=s,
                    endpoint_to_host_arcsec=dist,within_one_arcsec=bool(dist<=1),**extended_profile(science,valid,local),
                   centerline_distance_arcsec=cross,axis_difference_deg=angle,**row))
    for filepath in sorted(images.glob('legacy_rbh1*fits')):
        with fits.open(filepath) as hdul:
            h=hdul[0];w=WCS(h.header).celestial
            for i,band in enumerate(h.header['BANDS']):
                # Middle of published wake. Avoid host and tip. Absolute calibration unnecessary.
                axial,across=np.meshgrid(np.linspace(2,7,40),np.linspace(-.4,.4,7))
                def strip(shift):
                    sky=HOST.directional_offset_by(PA,axial*u.arcsec).directional_offset_by(
                            PA+90*u.deg,(shift+across)*u.arcsec)
                    x,y=w.world_to_pixel(sky)
                    return np.mean(map_coordinates(h.data[i],[y,x],order=1))
                def contrast(shift):
                    return strip(shift)-.5*(strip(shift-1.5)+strip(shift+1.5))
                parallel=np.array([contrast(s) for s in np.r_[-24:-4:1.25,5:25:1.25]])
                central=contrast(0); med=np.median(parallel);scatter=1.4826*np.median(abs(parallel-med))
                report['Legacy'].append(dict(release=h.header['VERSION'],band=band,contrast=float(central),
                    units='native image units/pixel',parallel_scatter=float(scatter),empirical_score=float((central-med)/scatter),
                    parallel_positions=len(parallel)))
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('images',type=Path);p.add_argument('results',type=Path);args=p.parse_args()
    (args.results/'rbh1_measured.json').write_text(json.dumps(measure(args.images,args.results),indent=2)+'\n')
