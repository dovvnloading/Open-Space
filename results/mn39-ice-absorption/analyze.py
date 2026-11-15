"""Replay the MN39 spectral and photometric measurements from offline inputs."""
import csv
import json
from datetime import date, timedelta
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
INPUT = ROOT/'inputs'
TARGET = (242.61055777065044, -51.357039296633985)
TARGET_ID = 15399993000757
# Band center and continuum anchors in observed micrometers. All fall within one detector.
BANDS = {'water': (3.033, (2.646, 2.794), (3.293, 3.476)),
         'carbon_dioxide': (4.257, (4.148, 4.184), (4.331, 4.368)),
         'carbon_monoxide': (4.669, (4.566, 4.635), (4.704, 4.773))}


def json_table(name):
    table = json.loads((INPUT/name).read_text())
    columns = [item['name'] for item in table['metadata']]
    return [dict(zip(columns, row)) for row in table['data']]


def write_table(name, rows):
    with (ROOT/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def spectrum(cube, wcs, ra, dec, radius, background):
    x, y = wcs.all_world2pix([[ra, dec]], 0)[0]
    yy, xx = np.indices(cube['image'].shape[1:]); dist = np.hypot(xx-x, yy-y)
    source = dist <= radius
    sky = (dist >= background[0]) & (dist <= background[1])
    values = []
    for image, hits in zip(cube['image'], cube['hits']):
        observed = source & np.isfinite(image) & (hits>0)
        annulus = sky & np.isfinite(image) & (hits>0)
        if observed.sum() != source.sum() or not annulus.any():
            values.append(np.nan)
        else:
            values.append(np.sum(image[observed] - np.median(image[annulus])))
    return np.asarray(values), source


def band_measurement(wavelength, signal, band):
    center, left, right = band
    centers, flux = [], []
    for lo, hi in (left, right):
        keep = (wavelength >= lo-.0005) & (wavelength <= hi+.0005)
        keep &= np.isfinite(signal)
        if np.count_nonzero(keep) < 2:
            return np.nan, np.nan
        centers.append(np.median(wavelength[keep])); flux.append(np.median(signal[keep]))
    idx = np.argmin(np.abs(wavelength-center))
    if min(flux) <= 0 or signal[idx] <= 0:
        return np.nan, np.nan
    continuum = np.exp(np.interp(wavelength[idx], centers, np.log(flux)))
    return signal[idx]/continuum, continuum


def main():
    cube = np.load(INPUT/'spherex.npz'); wave = cube['wavelength'].astype(float)
    wcs = WCS(fits.Header.fromstring(cube['header'].item()))
    standard, source = spectrum(cube, wcs, *TARGET, 1.5, (3, 5))
    records = []
    for i,w in enumerate(wave):
        f = standard[i]
        records.append(dict(channel=i+1, wavelength_um=round(w,3),
                            aperture_sum_mjy_sr='' if not np.isfinite(f) else round(float(f),5),
                            min_hits=int(cube['hits'][i][source].min()),
                            flags_or=int(np.bitwise_or.reduce(cube['flags'][i][source]))))
    ratios=[]
    for radius in (1.0,1.5,2.0):
        for annulus in ((3,5),(4,6)):
            spect,_ = spectrum(cube,wcs,*TARGET,radius,annulus)
            for label,band in BANDS.items():
                ratio,continuum = band_measurement(wave,spect,band)
                ratios.append(dict(feature=label,radius_pix=radius,background_pix=f'{annulus[0]}-{annulus[1]}',
                                   transmission=round(float(ratio),5),continuum_mjy_sr_sum=round(float(continuum),4)))
    # Compare nearby bright catalog sources selected without using spectral depths.
    candidates = ascii.read(INPUT/'allwise-field.tbl', format='ipac')
    comparisons=[]
    for row in candidates:
        if not (45 < row['dist'] < 110 and row['w1mpro'] < 10 and row['w2mpro'] < 10):
            continue
        s,field_mask = spectrum(cube,wcs,row['ra'],row['dec'],1.5,(3,5))
        entry={'source':row['designation'], 'distance_arcsec':round(float(row['dist']),2)}
        for label,band in BANDS.items():
            ratio,_ = band_measurement(wave,s,band)
            entry[label]='' if not np.isfinite(ratio) else round(float(ratio),4)
        comparisons.append(entry)
        if row['designation'] == 'J161036.07-512153.8':
            for radius in (1.0,1.5,2.0):
                for annulus in ((3,5),(4,6)):
                    neighbor,_=spectrum(cube,wcs,row['ra'],row['dec'],radius,annulus)
                    for label,band in BANDS.items():
                        ratio,continuum=band_measurement(wave,neighbor,band)
                        ratios.append(dict(feature=label,radius_pix=radius,background_pix=f'{annulus[0]}-{annulus[1]}',
                                           transmission=round(float(ratio),5),continuum_mjy_sr_sum=round(float(continuum),4),
                                           source='comparison'))
            for i,w in enumerate(wave):
                records[i]['neighbor_aperture_sum_mjy_sr']='' if not np.isfinite(s[i]) else round(float(s[i]),5)
                records[i]['neighbor_min_hits']=int(cube['hits'][i,field_mask].min())
                records[i]['neighbor_flags_or']=int(np.bitwise_or.reduce(cube['flags'][i,field_mask]))
    for row in ratios:
        row.setdefault('source','MN39')
    write_table('band-tests.csv',ratios)
    write_table('spectrum.csv',records)
    write_table('comparison-spectra.csv',comparisons)
    # Treat calendar years as epochs for VISTA, preserving the individual exposure counts.
    grouped={}
    for point in json_table('virac2-photometry.json'):
        if point['filter']!='Ks' or point['phot_flag'] or point['ambiguous_match'] or point['emag']>=.12:
            continue
        year=(date(1858,11,17)+timedelta(days=int(point['mjdobs']))).year
        grouped.setdefault((point['sourceid'],year),[]).append(point['mag'])
    vista=[dict(source_id=sid,year=year,count=len(vals),median_ks=round(float(np.median(vals)),4))
           for (sid,year),vals in sorted(grouped.items())]
    write_table('vista-epochs.csv',vista)
    table=ascii.read(INPUT/'neowise-target.tbl', format='ipac', converters={'moon_masked':str,'cc_flags':str})
    table.sort('mjd'); epochs=[]; groups=[]
    for row in table:
        if (row['cc_flags'][:2]!='00' or row['moon_masked']!='00' or row['dist']>=2 or row['qual_frame']<5
            or any(np.ma.is_masked(row[k]) or row[k]>=.12 for k in ('w1sigmpro','w2sigmpro'))):
            continue
        if not groups or row['mjd']-groups[-1][-1]['mjd']>40:
            groups.append([])
        groups[-1].append(row)
    for group in groups:
        if len(group)>=5:
            epochs.append(dict(n=len(group),**{k:round(float(np.median([r[k] for r in group])),4) for k in ('mjd','w1mpro','w2mpro')}))
    write_table('wise-epochs.csv',epochs)
    print('Aperture 1.5 px, annulus 3–5 px:',[r['transmission'] for r in ratios if r['radius_pix']==1.5 and r['background_pix']=='3-5'])

if __name__=='__main__':
    main()
