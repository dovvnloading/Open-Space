"""Circular counterfactual ephemerides and necessary bound-motion tests at 500–1000 au."""
from pathlib import Path
import argparse
import csv
import gzip
import json

import erfa
import numpy as np

HERE = Path(__file__).resolve().parent
MU = 0.0002959122082855911  # solar GM, au^3/day^2
C = 173.144632674240  # au/day
ARCSEC = 180 * 3600 / np.pi


def unit(ra, dec):
    return np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], axis=-1)


def angle(u, v):
    return np.arctan2(np.linalg.norm(np.cross(u, v), axis=-1), np.sum(u*v, axis=-1))


def tt(mjd):
    return erfa.taitt(*erfa.utctai(2400000.5, mjd))


def earth(mjd):
    return erfa.epv00(*tt(mjd))[0]["p"]


def read_astrometry():
    rows = {}
    with (HERE/'astrometry.csv').open() as f:
        for row in csv.DictReader(f):
            ra = 15*(float(row['ra_hours']) + float(row['ra_minutes'])/60 + float(row['ra_seconds'])/3600)
            dec = float(row['dec_degrees']) - float(row['dec_arcmin'])/60 - float(row['dec_arcsec'])/3600
            rows[row['source'], row['label']] = (float(row['utc_mjd']), *np.radians([ra,dec]))
    return rows


def circular_state(record, radii, headings):
    t, ra, dec = record
    n = unit(ra, dec)
    r, pa = np.broadcast_arrays(radii, np.radians(headings))
    e = earth(t)
    projection = e @ n
    delta = np.sqrt(projection**2 + r**2 - e @ e) - projection
    pos = e + delta[...,None]*n
    pole = pos/r[...,None]
    east = np.array([-np.sin(ra),np.cos(ra),0])
    east = east - (pole @ east)[...,None]*pole
    east /= np.linalg.norm(east,axis=-1)[...,None]
    north = np.cross(pole,east)
    tangent = np.cos(pa)[...,None]*north + np.sin(pa)[...,None]*east
    # Store state at emission. Planetary positions are heliocentric.
    return pos, tangent, r, np.sqrt(MU/r**3), sum(tt(t))-2400000.5-delta/C


def predict(state, mjd):
    pos, tangent, r, w, t0 = state
    e = earth(mjd)
    dt = sum(tt(mjd))-2400000.5-t0
    delay = r/C
    for _ in range(3):
        phi = w*(dt-delay)
        d = np.cos(phi)[...,None]*pos + (r*np.sin(phi))[...,None]*tangent - e
        delta = np.linalg.norm(d,axis=-1)
        delay = delta/C
    return d/delta[...,None]


def bound_test(start, stop, date_padding=0):
    t1,ra1,dec1 = start
    t2,ra2,dec2 = stop
    n1,n2=unit(ra1,dec1),unit(ra2,dec2)
    chord=earth(t2)-earth(t1)
    projected = np.linalg.norm(chord-(n1@chord)*n1)
    theta=angle(n1,n2)
    # Margin covers printed centroid differences/systematics. Gna gives only dates,
    # so allow independent times anywhere on those days (combined +/-1 day).
    projected_lower=projected-date_padding*0.018
    max_observed_angle=theta+2/ARCSEC
    required=projected_lower-1002*np.sin(max_observed_angle)
    # Delay difference bounded by (1002-498)/c < 3 d for the full distance bracket.
    # A bound orbit in this bracket stays far outside 495 au on these short spans.
    max_displacement=np.sqrt(2*MU/495)*(abs(t2-t1)+date_padding+3)
    return dict(interval_days=abs(t2-t1), observed_arcsec=theta*ARCSEC,
                projected_earth_chord_au=projected, extra_allowed_astrometric_error_arcsec=2,
                combined_date_padding_days=date_padding, required_displacement_lower_bound_au=required,
                bound_displacement_upper_bound_au=max_displacement)


def main():
    rows=read_astrometry()
    radii=np.array([500,700,1000])[:,None]
    headings=np.arange(0,360,0.5)[None,:]
    state=circular_state(rows['U','discovery_2016_fig5_table5'], radii, headings)
    grid_min={}
    for label in ['2014_band7','2014_band9','2016_band6','2017_band8','2017_band7']:
        t,ra,dec=rows['U',label]
        residual=angle(predict(state,t), unit(ra,dec))*ARCSEC
        grid_min[label]={str(r):round(float(val),3) for r,val in zip(radii.ravel(),np.min(residual,axis=1))}
    parser=argparse.ArgumentParser()
    parser.add_argument('--full',action='store_true',help='also write the four individual-track tables')
    args=parser.parse_args()
    with (HERE/'counterfactual_monthly.csv').open('w') as summary:
        writer=csv.writer(summary)
        writer.writerow(['utc_date','r_au','headings','ra_min','ra_max','dec_min','dec_max'])
        for year in range(2023,2027):
            gz = None
            if args.full:
                gz = gzip.GzipFile(filename=str(HERE/f'counterfactual_{year}.csv.gz'),mode='wb',mtime=0)
                gz.write(b'utc_date,r_au,heading_east_of_north_deg,ra_deg,dec_deg\n')
            try:
                for month in range(1,13):
                    mjd=erfa.cal2jd(year,month,1)[1]
                    d=predict(state,mjd)
                    ra=np.degrees(np.arctan2(d[...,1],d[...,0])) % 360
                    dec=np.degrees(np.arcsin(d[...,2]))
                    date=f'{year}-{month:02d}-01'
                    for ir,r in enumerate(radii.ravel()):
                        bounds=[ra[ir].min(), ra[ir].max(), dec[ir].min(), dec[ir].max()]
                        writer.writerow([date,r,headings.size,*np.round(bounds,7)])
                        if gz:
                            for ih,pa in enumerate(headings.ravel()):
                                gz.write(f'{date},{r},{pa:.1f},{ra[ir,ih]:.7f},{dec[ir,ih]:.7f}\n'.encode())
            finally:
                if gz:
                    gz.close()
    checks={
        'U_2016_Dec_2017_Mar': bound_test(rows['U','2016_band6'],rows['U','2017_band7']),
        'Gna_two_claimed_detections': bound_test(rows['Gna','discovery'],rows['Gna','second'],1),
    }
    report=dict(U_min_residual_arcsec=grid_min,necessary_bound_tests=checks,
                discovery_original_vs_reduced_arcsec=angle(unit(*rows['U','discovery_2016_fig5_table5'][1:]),unit(*rows['U','reanalysis_2015_band8'][1:]))*ARCSEC)
    (HERE/'result.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
