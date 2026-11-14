"""Register the four independent XMM pointings to the local Chandra frame."""
import csv
import json
import math
import statistics
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def rows(name):
    with (ROOT / name).open(encoding='utf-8', newline='') as source:
        return list(csv.DictReader(source))


def coordinates(row, ra='ra', dec='dec'):
    return float(row[ra]), float(row[dec])


OPTICAL = rows('optical.csv')
AUGUST = coordinates(OPTICAL[3], 'ra_deg', 'dec_deg')


def offset(coord):
    """Local angular offsets east and north, in arcseconds, from the H-alpha source."""
    return ((coord[0] - AUGUST[0]) * math.cos(math.radians(AUGUST[1])) * 3600,
            (coord[1] - AUGUST[1]) * 3600)


def median_point(points):
    return tuple(map(statistics.median, zip(*points)))


def rounded(value):
    return round(value, 3)


def register(xmm, chandra, radius):
    matches = []
    for entry, pos in xmm:
        distance = math.hypot(*pos)
        if not 12 < distance <= radius * 60:
            continue
        if float(entry['ep_extent']) or float(entry['ep_det_ml']) < 25 or float(entry['radec_err']) > 2:
            continue
        if any(math.dist(pos, other) < 9 for other_entry, other in xmm if entry != other_entry):
            continue
        nearest, second = sorted((math.dist(pos, ref), ref) for ref in chandra)[:2]
        if nearest[0] <= 5 and second[0] >= 8:
            matches.append(tuple(r - p for r, p in zip(nearest[1], pos)))
    for _ in range(3):
        center = median_point(matches)
        matches = [move for move in matches if math.dist(move, center) < 2.75]
    return len(matches), median_point(matches), tuple(map(statistics.pstdev, zip(*matches)))


def main():
    detections = rows('xmm-field.csv')
    chandra = [offset(coordinates(r)) for r in rows('chandra-field.csv')]
    neighbors = rows('neighbor.csv')
    positions = {r['observation']: offset(coordinates(r, 'ra_deg', 'dec_deg')) for r in OPTICAL}
    positions['PS1 neighbor'] = offset(coordinates(neighbors[0], 'raMean', 'decMean'))
    report = {'optical_may_to_august_days': rounded(float(OPTICAL[3]['mjd']) - float(OPTICAL[0]['mjd'])),
              'chandra_sources': len(chandra), 'xmm_field_detections': len(detections), 'pointings': []}
    # Sasaki et al. (2018), Table A.1, source 53. Their reduction of the same XMM field
    # used optical/infrared reference positions, independently of this Chandra registration.
    published = offset((15 * (43 / 60 + 47.788 / 3600), 41 + 33 / 60 + 12.38 / 3600))
    report['sasaki_53_separation_arcsec'] = {k: rounded(math.dist(published, v)) for k, v in positions.items()}
    neighbors = [r for r in rows('catalog-neighbors.csv') if r['name'] not in ('M31N 2014-08a', 'M31N 2015-06c')]
    separations = [math.hypot(*offset(coordinates(r, 'ra_deg', 'dec_deg'))) for r in neighbors]
    report['other_local_catalog_entries'] = {str(r): sum(d < r for d in separations) for r in (10, 30, 60, 120, 300)}
    for o in sorted(rows('observations.csv'), key=lambda item: item['start_utc']):
        oid = o['observation_id']
        xmm = [(r, offset(coordinates(r))) for r in detections if r['observation_id'] == oid]
        entry, raw = min(xmm, key=lambda pair: math.hypot(*pair[1]))
        mjd = (datetime.fromisoformat(o['start_utc']) - datetime(1858, 11, 17)).total_seconds() / 86400
        out = {'obs_id': oid, 'source': int(entry['src_num']), 'utc': o['start_utc'],
               'days_since_may_6': rounded(mjd - float(OPTICAL[0]['mjd'])),
               'ep_det_ml': rounded(float(entry['ep_det_ml'])),
               'hr1': rounded(float(entry['ep_hr1'])), 'hr2': rounded(float(entry['ep_hr2'])),
               'statistical_position_error_arcsec': rounded(float(entry['radec_err'])), 'registration': {}}
        for radius in (8, 12, 18):
            count, translation, spread = register(xmm, chandra, radius)
            corrected = tuple(p + t for p, t in zip(raw, translation))
            out['registration'][str(radius)] = {'references': count, 'shift_east_north_arcsec': list(map(rounded, translation)),
                'scatter_east_north_arcsec': list(map(rounded, spread)),
                'registered_east_north_arcsec': list(map(rounded, corrected)),
                'separation_arcsec': {k: rounded(math.dist(corrected, v)) for k, v in positions.items()}}
        report['pointings'].append(out)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
