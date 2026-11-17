"""Replay signed VLASS forced photometry and equal-radius offset trials offline."""
import csv
import hashlib
import json
import math
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent


def photometry(packet, at, controls):
    with fits.open(ROOT/packet['path']) as hdus:
        data, hdr = hdus[0].data, hdus[0].header
        if hashlib.sha256(data.tobytes()).hexdigest() != packet['data_sha256']:
            raise ValueError(f"Pixel checksum mismatch for {packet['path']}")
        w = WCS(hdr).celestial
        x, y = w.world_to_pixel(at)
        scale = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600
        yy, xx = np.indices(data.shape)
        rr = np.hypot(xx-x, yy-y) * scale
        annulus = data[(rr >= 12) & (rr <= 40) & np.isfinite(data)]
        beam_pixels = math.pi * hdr['BMAJ'] * hdr['BMIN'] * 3600**2 / (4 * np.log(2)*scale**2)
        if len(annulus) < 50*beam_pixels:
            raise ValueError(f"Fewer than 50 finite noise beams in {packet['path']}")
        med = np.median(annulus)
        rms = 1.4826022185 * np.median(abs(annulus-med))
        if not np.isfinite(rms) or rms <= 0:
            raise ValueError('Invalid local noise')
        points = [(-1,at)]
        if controls:
            points.extend((pa,at.directional_offset_by(pa*u.deg,25*u.arcsec)) for pa in range(0,360,45))
        fluxes = []
        for pa, point in points:
            px, py = np.rint(w.world_to_pixel(point)).astype(int)
            if 0 <= py < data.shape[0] and 0 <= px < data.shape[1] and np.isfinite(data[py,px]):
                flux = (data[py,px]-med)*1000
                fluxes.append((pa, flux, flux/(rms*1000)))
        return dict(rms_mjy_beam=rms*1000, median_mjy_beam=med*1000,
                    annulus_beams=len(annulus)/beam_pixels, fluxes=fluxes)


def table(name, rows, columns):
    with (ROOT/name).open('w', newline='') as out:
        w = csv.DictWriter(out, columns, lineterminator='\n'); w.writeheader(); w.writerows(rows)


def main():
    run = json.loads((ROOT/'extraction.json').read_text())
    targets, offsets, distinct = [], [], []
    for field in run['fields']:
        positions = []
        for group in field['groups']:
            for s in group['rows']:
                base = dict(field=field['field'], utc=field['obs_date'], sample=s['sample'], target=s['target'],
                            gaia_dr3_id=s['gaia_dr3_id'], ra_deg=s['ra'], dec_deg=s['dec'], pm_applied=s['pm_applied'])
                at = SkyCoord(s['ra'],s['dec'],unit='deg')
                v = photometry(group['V'],at,True); i = photometry(group['I'],at,False)
                vf = next((f for pa,*f in v['fluxes'] if pa == -1),[None,None]); inf = next((f for pa,*f in i['fluxes'] if pa == -1),[None,None])
                # CNS5 gives an alternative reference epoch for the same T dwarf (see README).
                alternate = s['sample'] == 'cns5' and s['target'] == 'CNS5 277'
                independently_counted = bool(not alternate and s['pm_applied'] and vf[0] is not None and
                    all(at.separation(prev).arcsec > field['beam_maj_arcsec'] for prev in positions))
                if independently_counted:
                    positions.append(at)
                on = dict(**base, independently_counted=independently_counted, alternate_astrometry=alternate,
                          v_mjy_beam=vf[0], v_mad_rms_mjy_beam=v['rms_mjy_beam'], v_snr=vf[1],
                          v_5sigma_mjy_beam=5*v['rms_mjy_beam'], v_background_mjy_beam=v['median_mjy_beam'],
                          v_annulus_beams=v['annulus_beams'], i_mjy_beam=inf[0], i_mad_rms_mjy_beam=i['rms_mjy_beam'], i_snr=inf[1])
                targets.append(on)
                if independently_counted:
                    distinct.append(on)
                offsets.extend(dict(**base, independently_counted=independently_counted, offset_pa_deg=pa, v_mjy_beam=flux, v_snr=snr)
                               for pa,flux,snr in v['fluxes'] if pa >= 0)
    if len(targets) == 0:
        raise RuntimeError('No photometry to replay')
    table('measurements.csv',targets,list(targets[0])); table('offsets.csv',offsets,list(offsets[0]))
    ctrl = [r for r in offsets if r['independently_counted']]
    sig = lambda rows: sum(bool(abs(r['v_snr']) >= 5) for r in rows)
    summary = dict(input_rows=run['samples'], files=len(run['fields']),
        label_counts={s:sum(r['sample']==s for r in targets) for s in run['samples']},
        pm_missing_label_rows=sum(not r['pm_applied'] for r in targets), independent_on=len(distinct),
        independent_controls=len(ctrl), n_5sigma_on=sig(distinct),n_5sigma_controls=sig(ctrl),
        max_abs_snr_on=max(abs(r['v_snr']) for r in distinct),max_abs_snr_offset=max(abs(r['v_snr']) for r in ctrl),
        five_sigma_mjy_beam_range=[min(r['v_5sigma_mjy_beam'] for r in distinct),max(r['v_5sigma_mjy_beam'] for r in distinct)])
    (ROOT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
