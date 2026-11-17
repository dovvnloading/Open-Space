"""Crossmatch native NEOWISE detections and static catalogues for each W2 residual."""
import csv
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time

HERE = Path(__file__).resolve().parent


def coordinates(rows):
    return SkyCoord([float(row['ra']) for row in rows], [float(row['dec']) for row in rows], unit='deg')


def catalog_near(label, visit, coord):
    with (HERE/f'{label}{visit}.response').open() as fin:
        rows = list(csv.DictReader(fin))
    sep = coord.separation(coordinates(rows)).arcsec
    return [dict(rows[i], separation_arcsec=round(sep[i], 3)) for i in np.argsort(sep) if sep[i] < 35]


def main():
    positions = json.loads((HERE/'positions.json').read_text())
    individual = list(csv.DictReader((HERE/'individual.csv').open()))
    result = []
    for track in positions:
        group = []
        for visit in track['visits']:
            v = visit['visit']
            pointing = coordinates([visit])[0]
            single = Table.read(HERE/f'l1visit{v}.response', format='ascii.ipac')
            native = [r for r in individual if int(r['template']) == track['template'] and int(r['visit']) == v]
            closest = []
            for row in native:
                _, scan, frame = row['frame'].split('_')
                subset = single[(single['scan_id'] == scan) & (single['frame_num'] == int(frame))]
                item = dict(frame=row['frame'], number_catalogued_in_cone=len(subset))
                if len(subset):
                    sep = coordinates([row])[0].separation(coordinates(subset)).arcsec
                    clean = (subset['w2snr'].filled(0) if hasattr(subset['w2snr'],'filled') else subset['w2snr']) >= 5
                    clean &= np.array([str(x)[1]=='0' for x in subset['cc_flags']])
                    i = np.argmin(sep)
                    item.update(nearest_arcsec=round(sep[i],3), nearest_W2_snr=float(subset['w2snr'][i]),
                                nearest_flags=str(subset['cc_flags'][i]),
                                clean_SNR5_within_3arcsec=int(np.sum((sep<=3)&clean)))
                closest.append(item)
            jpl = json.loads((HERE/f'jpl{v}.response').read_text())
            wise = catalog_near('allwise',v,pointing)
            for item in wise:
                gaia = catalog_near('gaia',v,coordinates([item])[0])
                if gaia and gaia[0]['separation_arcsec'] <= 1:
                    item['gaia_nearest_at_catalog_epoch'] = gaia[0]
            group.append(dict(visit=v, utc=Time(visit['mjd'],format='mjd').isot,
                              native_valid_science=sum(np.isfinite(float(r['science_filtered'])) for r in native).item(),
                              allwise_within_35arcsec=wise,
                              gaia_dr3_at_catalog_epoch_within_35arcsec=catalog_near('gaia',v,pointing),
                              selected_l1b_frames=closest,jpl_known_object_outcome=jpl.get('warning',jpl.get('n_second_pass'))))
        coords=coordinates(track['visits'])
        result.append(dict(template=track['template'], visits=group,
                           predicted_steps_arcmin=[float(coords[i].separation(coords[j]).arcmin) for i,j in [(0,1),(1,2),(0,2)]]))
    (HERE/'vetting.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
