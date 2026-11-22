"""Blind pixel selection and calibration on each complete, checksum-verified release."""
import argparse
from collections import Counter, defaultdict
import gzip
import json
import multiprocessing as mp
import os
from pathlib import Path
import random
import time

import numpy as np

from acquisition import Spectra, prepare
import spectral as spec

# Vacuum wavelengths including common stellar/chromospheric, nebular and UV features.
# An unmatched peak at z_cat is a triage result. A different redshift or a blend can fit it.
REST = np.sort(np.array([
    912, 972.537, 977.02, 1025.722, 1031.926, 1037.617, 1084.6, 1175.7, 1190.4,
    1193.3, 1200, 1206.50, 1215.6701, 1238.821, 1242.804, 1260.42, 1264.7,
    1302.17, 1304.37, 1309.28, 1334.53, 1335.71, 1393.755, 1397.2, 1402.770,
    1404.8, 1483.32, 1486.50, 1526.71, 1533.43, 1548.204, 1550.781, 1608.45,
    1640.42, 1660.809, 1666.150, 1670.79, 1718, 1750, 1854.72, 1862.79, 1882.71,
    1892.03, 1906.683, 1908.734, 2143, 2326, 2424, 2471, 2586.65, 2600.17,
    2626.45, 2796.35, 2803.53, 2852.96, 3133, 3188, 3204, 3346.79, 3426.85,
    3445, 3587, 3615, 3635, 3672, 3683, 3692, 3698, 3705, 3713, 3722,
    3727.092, 3729.875, 3735, 3751.2, 3771.7, 3799, 3820, 3836.5, 3869.86,
    3889.75, 3890.2, 3934.78, 3968.59, 3969.59, 3971.19, 4009, 4027.3,
    4069.75, 4077.5, 4102.892, 4145, 4201, 4228, 4268, 4341.684, 4364.44,
    4389, 4472.74, 4543, 4641, 4651, 4659.38, 4687.02, 4712.6, 4741.5,
    4862.683, 4923.3, 4960.295, 4981.7, 5008.240, 5017, 5049, 5159, 5170, 5199,
    5202, 5272, 5413, 5519, 5539, 5877, 6302, 6366, 6549.86, 6564.61, 6585.27]))
CAT_KEYS = ('release','source_id','detectid','shotid','ra','dec','source_type','z_hetdex',
            'z_hetdex_src','z_hetdex_conf','det_type','wave','field','ifuslot','flag_best',
            'agn_flag','sigma','ra_det','dec_det','cls_diagnose','plya_classification','fwhm','sn','p_conf','p_cnn','major','minor','gmag')


def json_default(obj):
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(str(type(obj)))


def dump(row):
    return json.dumps(row, default=json_default, allow_nan=False, separators=(',', ':'))


def jsonlines(path):
    with gzip.open(path, 'rt') as src:
        for line in src:
            yield json.loads(line)


def channels(wave, flux, error, peak, z):
    """Apply identical hypotheses to source rows and their shuffled counterparts."""
    if peak['single_pixel'] or peak['sky_near'] or peak.get('badcal', False):
        return []
    found = []
    if peak['ew_extreme'] and 1.88 <= peak['z_lya'] <= 3.52:
        found.append('ew')
    zh = peak['z_heii']
    if zh is not None:
        anchored = z is not None and abs(z - zh) < 0.008
        if anchored or (peak.get('lya_pre_snr') or 0) >= 4.0:
            found.append('heii')
    if z is not None and z >= 0 and peak['snr'] >= 9:
        red = 1 + z
        tolerance = max(8, 2*peak['sigma_A'], peak['wave_vac']*600/spec.C)
        if np.min(np.abs(REST*red - peak['wave_vac'])) > tolerance:
            found.append('unknown')
    return found


def screen_chunk(wave, flux, error, redshifts, offset, kind):
    counts = Counter(spectra=len(flux), usable_spectra=int(np.sum(np.sum(np.isfinite(error) & (error>0), axis=1) >= 32)))
    rows = []
    for index, peak in spec.search_batch(wave, flux, error, density_unit=1e-17):
        counts['peaks'] += 1
        if peak['single_pixel']:
            counts['single_pixel_flags'] += 1
        if peak['sky_near'] or peak.get('badcal', False):
            counts['sky_or_calibration_flags'] += 1
        z = float(redshifts[index])
        kinds = channels(wave, flux[index], error[index], peak, z if np.isfinite(z) and z > -0.01 else None)
        for channel in kinds:
            counts[channel] += 1
            rows.append(dict(peak, channel=channel, row_index=index+offset, stream=kind))
    return counts, rows


CAT = None
BATCH = 256


def work(start):
    cat = CAT
    flux, error = cat.read(start, min(start+BATCH, cat.size))
    z = cat.info['z_hetdex'][start:start+len(flux)]
    counts, rows = screen_chunk(cat.wave, flux, error, z, start, 'observed')
    control = Counter()
    fake_rows = []
    # One randomly assigned block per 8-block stratum. No decision uses its contents.
    block_no = start//BATCH
    if len(flux) > 1 and block_no % 8 == random.Random(700001 + block_no//8).randrange(8):
        fy, fe = spec.shuffle_fibres(flux, error, seed=410087 + start)
        control, fake_rows = screen_chunk(cat.wave, fy, fe, z, start, 'shuffled')
    injections = []
    # One random row per eight blocks, six independent trials per selected spectrum.
    if block_no % 8 == random.Random(230009 + block_no//8).randrange(8):
        i = random.Random(930829 + start).randrange(len(flux))
        trials = spec.simulate(cat.wave, flux[i], error[i], seed=502013+start,
                               n_inject=6, density_unit=1e-17)['injections']
        injections = [dict(trial, row_index=start+i) for trial in trials]
    return counts, rows, control, fake_rows, injections


def scan(cat, output, workers):
    global CAT
    CAT = cat
    observed, controls, inj_counts = Counter(), Counter(), Counter()
    stats = {}
    step = time.monotonic()
    with gzip.open(output/'screen.jsonl.gz','wt',compresslevel=3) as obs, \
            gzip.open(output/'shuffled.jsonl.gz','wt',compresslevel=3) as ctl, \
            gzip.open(output/'injections.jsonl.gz','wt',compresslevel=3) as inj, \
            mp.get_context('fork').Pool(workers) as pool:
        results = pool.imap(work, range(0,cat.size,BATCH), chunksize=2)
        for counts, rows, null_counts, fake_rows, injections in results:
            observed.update(counts)
            controls.update(null_counts)
            for fh, data in [(obs,rows),(ctl,fake_rows),(inj,injections)]:
                for row in data:
                    fh.write(dump(row)+'\n')
            for row in injections:
                inj_counts[row['status']] += 1
                if row['status']=='injected' and not row['was_blended']:
                    inj_counts['unblended'] += 1
                    inj_counts['unblended_recovered'] += int(row['recovered'])
            if time.monotonic()-step>60:
                print(cat.release, observed, flush=True)
                step=time.monotonic()
    stats.update(release=cat.release, observed=dict(observed), shuffled=dict(controls),
                 injections=dict(inj_counts), distinct_release_source_ids=int(len(np.unique(cat.src_ids))),
                 input_rows=cat.size, wavelength_air_A=[float(cat.wave[0]),float(cat.wave[-1])])
    (output/'scan-summary.json').write_text(json.dumps(stats,indent=2)+'\n')


def followup(cat, out, cache):
    from vetting import spatial, vet, blank_apertures
    from vet_cubes import halo_summary
    from positive_controls import positive_controls
    from ew_rescue import rescue
    from calibration import summarize

    positive_controls(cat,out)
    summarize(out)
    blanks=Counter()
    with gzip.open(out/'offsource.jsonl.gz','wt') as fh:
        for meta, flux, err in blank_apertures(cat,cache,n_cubes=8):
            counts, rows = screen_chunk(cat.wave,flux,err,np.full(len(flux),np.nan),0,'offsource') if len(flux) else (Counter(), [])
            blanks.update(counts)
            fh.write(dump(dict(input=meta,counts=counts,features=rows))+'\n')
    halo_summary(out)
    (out/'offsource-summary.json').write_text(json.dumps(blanks,indent=2)+'\n')
    from itertools import chain
    from repeats import repeat_search
    from pair_search import pair_search
    paired = pair_search(cat,out)
    repeats = repeat_search(cat,out)
    imaging = rescue(cat,out)
    groups = spatial(cat)
    print('Spatial prescreen',len(groups),'repeated-visit leads',len(repeats),
          'imaging-led objects',len(imaging),flush=True)
    # vet() retrieves complete metadata by row after blind pixel selection.
    vet(chain(jsonlines(out/'screen.jsonl.gz'),repeats,imaging,groups,paired),cat,out)
    from control_followup import control_followup
    control_followup(cat,out)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('release',choices=['sc1','sc2'])
    ap.add_argument('mode',choices=['scan','vet'])
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=3)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    path=prepare(args.release,args.cache,args.output)
    with Spectra(path,args.release) as cat:
        if args.mode=='scan':
            scan(cat,args.output,max(1,args.workers))
        else:
            followup(cat,args.output,args.cache)


if __name__=='__main__':
    main()
