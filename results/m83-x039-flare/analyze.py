"""Reproduce the M83 X039 measurements from the public event file and catalogs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.special import betainc
from scipy.stats import rankdata, spearmanr

RA, DEC = 204.190875, -29.830555556
EVENT_URL = 'https://cxc.cfa.harvard.edu/cdaftp/byobsid/2/13202/primary/acisf13202N003_evt2.fits.gz'
SHA256 = 'ab40ebbc7071a3a921f4c3a47f92b59ed1c095920a01df1e6dd37263addceda8'
NSC_ID = '147490_8194'
PS_ID = '72202041909333633'


def catalog_urls():
    nsc = 'https://datalab.noirlab.edu/scs/nsc_dr2/object?'
    ps = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?'
    return {
        'chandra.fits.gz': EVENT_URL,
        'nsc.xml': nsc + urlencode(dict(RA=RA, DEC=DEC, SR=65 / 3600, VERB=3)),
        'panstarrs.csv': ps + urlencode(dict(ra=RA, dec=DEC, radius=60 / 3600, pagesize=12000,
            columns='[objID,raMean,decMean,nDetections,zMeanPSFMag,yMeanPSFMag]')),
    }


def fetch_inputs(cache):
    cache.mkdir(parents=True, exist_ok=True)
    for name, url in catalog_urls().items():
        path = cache / name
        if not path.exists():
            with urlopen(url, timeout=180) as response:
                body = response.read()
            temporary = path.with_suffix('.part')
            temporary.write_bytes(body)
            temporary.replace(path)
    if hashlib.sha256((cache / 'chandra.fits.gz').read_bytes()).hexdigest() != SHA256:
        raise ValueError('Unexpected Chandra file digest')


def sky_wcs(header):
    wcs = WCS(naxis=2)
    for field, key in [('ctype', 'TCTYP'), ('crval', 'TCRVL'), ('cdelt', 'TCDLT'), ('crpix', 'TCRPX')]:
        setattr(wcs.wcs, field, [header[key + '11'], header[key + '12']])
    return wcs


def angular_distance(ra, dec):
    # Tangent-plane approximation over at most 65 arcsec.
    return 3600 * math.hypot((float(ra) - RA) * math.cos(math.radians(DEC)), float(dec) - DEC)


def photons(cache, output):
    with fits.open(cache / 'chandra.fits.gz') as hdus:
        events = hdus['EVENTS']
        header = events.header
        xy = sky_wcs(header).all_world2pix([[RA, DEC]], 1)[0]
        off_axis = np.linalg.norm(xy - [header['TCRPX11'], header['TCRPX12']]) * 0.492 / 60
        radius = max(1.75, 1.3 + 0.07 * off_axis**2)
        table = events.data
        distance = np.hypot(table['x'] - xy[0], table['y'] - xy[1]) * 0.492
        chip = 7
        gti = next(hdu.data for hdu in hdus if hdu.name == 'GTI' and hdu.header['CCD_ID'] == chip)
        if len(gti) != 1:
            raise ValueError('The analysis assumes a single continuous good-time interval')
        start, stop = float(gti['START'][0]), float(gti['STOP'][0])
        duration = stop - start
        good = ((table['ccd_id'] == chip) & (table['energy'] >= 350) & (table['energy'] <= 8000)
                & (table['time'] >= start) & (table['time'] < stop))
        source = table[good & (distance < radius)]
        source = source[np.argsort(source['time'])]
        time = source['time'].astype(float) - start
        energy = source['energy'].astype(float)
        annulus = good & (distance > 5) & (distance <= 12)
        background_time = table['time'][annulus] - start
        background_energy = table['energy'][annulus]
        n = len(time)
        with (output / 'photons.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['seconds_from_chip7_gti_start', 'energy_eV', 'chipx', 'chipy'])
            writer.writerows((f'{t:.6f}', f'{e:.3f}', int(row['chipx']), int(row['chipy']))
                             for t, e, row in zip(time, energy, source))
        # The screen considered every contiguous run at each allowed photon count.
        m = 32
        lag = time[m - 1:] - time[:n - m + 1]
        j = int(np.argmin(lag))
        fraction = (lag[j] + 3.241) / duration
        p_window = float(betainc(m - 1, n - m + 2, fraction))
        p_scan = min(1, p_window * (n - m + 1) * 3943 * 12)
        x = rankdata(time)
        x -= x.mean()
        y = rankdata(energy)
        y -= y.mean()
        observed = abs(x @ y)
        generator = np.random.default_rng(13202)
        more_extreme = sum(abs(x @ generator.permutation(y)) >= observed for _ in range(100000))
        windows = [(0, 55000), (55000, 65000), (65000, 79000), (79000, duration)]
        blocks = []
        for lo, hi in windows:
            signal = (time >= lo) & (time < hi)
            local = (background_time >= lo) & (background_time < hi)
            blocks.append(dict(start_s=lo, stop_s=round(hi, 3), source=int(signal.sum()),
                annulus=int(local.sum()), expected_background=round(local.sum()*radius**2 / (12**2-5**2), 3),
                source_median_eV=round(float(np.median(energy[signal])), 1)))
        reduced = {}
        for test_radius in [1, 1.5, 1.75, 2, 2.5, 3]:
            mask = good & (distance < test_radius)
            reduced[str(test_radius)] = [int(mask.sum()), round(float(spearmanr(table['time'][mask], table['energy'][mask]).statistic), 4)]
        utc = Time(header['MJDREF'] + (start+time[j])/86400, format='mjd', scale='tt').utc.isot
        return dict(obsid=13202, file_sha256=SHA256, gti_start_tt_mission_s=start,
            gti_seconds=round(duration, 3), off_axis_arcmin=round(off_axis, 4), aperture_arcsec=round(radius, 4),
            count=n, annulus_5_12_arcsec_count=int(annulus.sum()), windows=blocks,
            shortest_32_photon_interval_s=round(float(lag[j]), 3), interval_start_utc=utc,
            conditional_constant_rate_p_family_bound=p_scan, time_energy_spearman=float(spearmanr(time, energy).statistic),
            time_energy_permutation_p_two_sided=(1+more_extreme)/100001,
            local_annulus_time_energy_spearman=float(spearmanr(background_time, background_energy).statistic),
            aperture_sensitivity_count_spearman=reduced)


def optical(cache):
    root = ElementTree.parse(cache / 'nsc.xml')
    fields = [field.attrib['name'] for field in root.findall('.//{*}FIELD')]
    rows = [dict(zip(fields, (cell.text for cell in row))) for row in root.findall('.//{*}TR')]
    target = next(row for row in rows if row['id'] == NSC_ID)
    controls = [row for row in rows if row['id'] != NSC_ID and int(row['ndet']) >= 30
                and 19.5 < float(row['imag']) < 21.5 and float(row['pmraerr']) < 5 and float(row['pmdecerr']) < 5]
    motion = np.array([[float(row[k]) for k in ('pmra', 'pmdec')] for row in controls])
    median = np.median(motion, axis=0)
    target_motion = np.array([float(target[k]) for k in ('pmra', 'pmdec')])
    rank = int(np.count_nonzero(np.linalg.norm(motion-median, axis=1) >= np.linalg.norm(target_motion-median)))
    with (cache / 'panstarrs.csv').open() as stream:
        ps = list(csv.DictReader(stream))
    counterpart = next(row for row in ps if row['objID'] == PS_ID)
    field = [row for row in ps if 10 < angular_distance(row['raMean'], row['decMean']) <= 60 and int(row['nDetections']) >= 2]
    separation = angular_distance(counterpart['raMean'], counterpart['decMean'])
    area = math.pi*(60**2-10**2)
    return dict(nsc_id=NSC_ID, nsc_offset_arcsec=round(angular_distance(target['ra'],target['dec']), 4),
        nsc={k: target[k] for k in ('ra','dec','ndet','rmag','rerr','imag','ierr','zmag','zerr','pmra','pmdec','pmraerr','pmdecerr','deltamjd')},
        nsc_motion_controls=len(controls), controls_at_least_as_large_motion=rank, local_median_motion=median.round(3).tolist(),
        ps1=counterpart, ps1_offset_arcsec=round(separation,4), ps1_control_count=len(field),
        local_poisson_p_within_offset=1-math.exp(-len(field)/area*math.pi*separation**2),
        local_poisson_p_within_csc_semimajor=1-math.exp(-len(field)/area*math.pi*0.803**2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('x039-cache'))
    parser.add_argument('--output', type=Path, default=Path('.'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    fetch_inputs(args.cache)
    result = {'xray': photons(args.cache, args.output), 'optical': optical(args.cache)}
    (args.output/'measurements.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
