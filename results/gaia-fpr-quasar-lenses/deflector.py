"""Optional catalogue deblends at independently fitted lens centres.

Use only on preselected systems. Legacy Surveys Tractor and PS1 are separate
catalogue reductions of the imaging; they are not separate observations from
the corresponding survey pixels. No deblend or morphology confirms lensing.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import time
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import requests
from astropy.table import Table
from scipy.optimize import linear_sum_assignment
from astropy.coordinates import angular_separation

CONE = 6.0  # arcseconds
CENTRE_RADIUS = 0.5
IMAGE_RADIUS = 0.4
PS1 = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.csv'
TAP = 'https://datalab.noirlab.edu/tap/sync'
PS1_COLUMNS = ['objID', 'raStack', 'decStack', 'raStackErr', 'decStackErr', 'raMean',
               'decMean', 'raMeanErr', 'decMeanErr', 'qualityFlag', 'primaryDetection',
               'bestDetection', 'nStackDetections']
for band in 'riz':
    PS1_COLUMNS += [band + s for s in ('PSFMag', 'PSFMagErr', 'KronMag',
                     'KronMagErr', 'psfQfPerfect', 'nFrames', 'psfMajorFWHM')]


def numeric(value):
    try:
        n = float(value)
        return n if math.isfinite(n) and n != -999 else None
    except (TypeError, ValueError):
        return None


def inverse_error(value, scale=1):
    return scale / math.sqrt(value) if value is not None and value > 0 else None


def separation(a, b):
    """Spherical separation from (ra,dec) to (ra,dec), in arcseconds."""
    return float(np.rad2deg(angular_separation(*np.deg2rad([*a, *b]))) * 3600)


def parse_tap(content):
    root = ElementTree.fromstring(content)
    statuses = [node.get('value', '').upper() for node in root.iter()
                if node.tag.endswith('INFO') and node.get('name', '').upper() == 'QUERY_STATUS']
    if not statuses or any(status != 'OK' for status in statuses):
        raise ValueError(f'TAP status {statuses}')
    table = Table.read(io.BytesIO(content), format='votable', use_names_over_ids=True)
    names = {k.lower(): k for k in table.colnames}
    if not {'ra', 'dec', 'type', 'brickname', 'objid'} <= set(names) or len(table) >= 201:
        raise ValueError('Incomplete or truncated Tractor cone')
    # Extract only published fields used below. Preserve the complete raw cone.
    selected = {'ra','dec','type','brickname','objid','release','brick_primary','ls_id',
                'ra_ivar','dec_ivar','shape_r','shape_r_ivar','maskbits','fitbits','ref_cat'}
    selected |= {key+'_'+band for key in ('nobs','flux','flux_ivar','psfsize','fracflux',
                                          'rchisq','fracmasked','fracin') for band in 'griz'}
    def cell(value):
        return value.decode().strip() if isinstance(value, bytes) else str(value).strip()
    return [{k: cell(row[names[k]]) for k in selected if k in names} for row in table]


def parse_ps1(content):
    table = csv.DictReader(io.StringIO(content.decode('utf-8-sig')))
    if not table.fieldnames or {s.lower() for s in PS1_COLUMNS} - {s.lower() for s in table.fieldnames}:
        raise ValueError('Malformed PS1 StackObjectView response')
    records = [{k.lower(): v for k, v in row.items()} for row in table]
    if len(records) >= 200:
        raise ValueError('PS1 result reached requested page size')
    return records


def get_cone(survey, centre, out):
    ra, dec = centre
    if survey == 'PS1_DR2_stack':
        url, parse = PS1, parse_ps1
        params = dict(ra=f'{ra:.10f}', dec=f'{dec:.10f}', radius=CONE/3600,
                      columns='['+','.join(PS1_COLUMNS)+']', pagesize=200)
    else:
        url, parse = TAP, parse_tap
        table = {'Legacy_DR10':'ls_dr10.tractor', 'Legacy_DR9':'ls_dr9.tractor'}[survey]
        query = (f"SELECT TOP 201 * FROM {table} AS g WHERE 1=CONTAINS("
                 f"POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{ra:.10f},{dec:.10f},{CONE/3600:.12f}))")
        params = dict(REQUEST='doQuery', LANG='ADQL', FORMAT='votable', QUERY=query, MAXREC=201)
    receipt = dict(url=url, parameters=params, radius_arcsec=CONE, attempts=[])
    for attempt in range(2):
        try:
            response = requests.get(url, params=params, timeout=(10, 35))
            receipt['attempts'].append(dict(status_code=response.status_code))
            response.raise_for_status()
            data = parse(response.content)
            path = out/(survey + ('.csv.gz' if survey == 'PS1_DR2_stack' else '.vot.gz'))
            path.write_bytes(gzip.compress(response.content, mtime=0))
            receipt.update(sha256=hashlib.sha256(response.content).hexdigest(),
                           bytes=len(response.content), raw_file=path.name, returned_rows=len(data))
            return data, receipt
        except Exception as exc:
            receipt['attempts'].append(dict(error=f'{type(exc).__name__}: {str(exc)[:300]}'))
            if attempt == 0:
                time.sleep(1)
    return None, receipt


def tractor(row):
    ra, dec = numeric(row['ra']), numeric(row['dec'])
    primary = row.get('brick_primary','true').lower()
    if None in (ra, dec) or primary in ('false', 'f', '0'):
        return None
    kind = row['type'].strip().upper()
    # A Tractor DUP is an explicit zero optical flux placeholder for Gaia.
    if kind == 'DUP':
        return None
    extended = kind in ('REX','EXP','DEV','SER')
    size = numeric(row.get('shape_r'))
    esize = inverse_error(numeric(row.get('shape_r_ivar')))
    bands = {}
    for band in 'griz':
        def get(key):
            return numeric(row.get(key+'_'+band))
        flux, ivar = get('flux'), get('flux_ivar')
        bands[band] = dict(nobs=get('nobs'), psf_fwhm_arcsec=get('psfsize'),
                           flux_nanomaggies=flux, flux_error_nanomaggies=inverse_error(ivar),
                           fracflux=get('fracflux'), fracmasked=get('fracmasked'), fracin=get('fracin'),
                           reduced_chi2=get('rchisq'))
    good_flux = any((b['nobs'] or 0)>0 and (b['flux_nanomaggies'] or 0)>0 and
                    b['flux_error_nanomaggies'] is not None and
                    b['flux_nanomaggies']/b['flux_error_nanomaggies'] >= 5 for b in bands.values())
    support = extended and good_flux and size is not None and esize is not None and size > 2*esize
    ident = row.get('ls_id') or '/'.join(row.get(k,'') for k in ('release','brickname','objid'))
    return dict(object_id=ident,ra=ra,dec=dec,type=kind,looks_extended=extended,
                morphology_quality=support,half_light_radius_arcsec=size,half_light_radius_error_arcsec=esize,
                ra_coord_error_arcsec=inverse_error(numeric(row.get('ra_ivar')),3600),
                dec_coord_error_arcsec=inverse_error(numeric(row.get('dec_ivar')),3600),
                maskbits=row.get('maskbits'),fitbits=row.get('fitbits'),ref_cat=row.get('ref_cat'),
                release=row.get('release'),bands=bands)


def panstarrs(row):
    def get(k):
        return numeric(row.get(k.lower()))
    ra, dec = get('raStack'), get('decStack')
    if None in (ra, dec):
        return None
    bands, good, extended = {}, [], []
    for band in 'riz':
        psf, kron = get(band+'PSFMag'), get(band+'KronMag')
        ep, ek = get(band+'PSFMagErr'), get(band+'KronMagErr')
        qf, frames = get(band+'psfQfPerfect'), get(band+'nFrames')
        delta = psf-kron if None not in (psf,kron) and 0 < min(psf,kron) < 40 and max(psf,kron) < 40 else None
        e_delta = math.hypot(ep,ek) if None not in (ep,ek) and min(ep,ek)>0 else None
        # PSF and Kron use the same pixels: quadrature is a nominal diagnostic;
        # their covariance is not supplied. Retain the separate errors.
        clean = (delta is not None and e_delta is not None and max(ep,ek)<=0.2 and
                 qf is not None and qf>=0.85 and frames is not None and frames>=2)
        if clean:
            good.append(band)
            if delta > max(0.1,3*e_delta):
                extended.append(band)
        bands[band] = dict(psf_mag=psf,kron_mag=kron,psf_mag_err=ep,kron_mag_err=ek,
                           psf_minus_kron_mag=delta,nominal_error_mag=e_delta,unmasked_fraction=qf,
                           frames=frames,psf_major_fwhm_arcsec=get(band+'psfMajorFWHM'))
    return dict(object_id=row['objid'],ra=ra,dec=dec,type='stack',looks_extended=bool(extended),
                morphology_quality=len(extended)>=2,valid_bands=good,extended_bands=extended,
                primary_detection=get('primaryDetection'),best_detection=get('bestDetection'),
                quality_flag=get('qualityFlag'),n_stack_detections=get('nStackDetections'),
                stack_ra_scatter_arcsec=get('raStackErr'),stack_dec_scatter_arcsec=get('decStackErr'),
                mean_ra_error_arcsec=get('raMeanErr'),mean_dec_error_arcsec=get('decMeanErr'),
                bands=bands)


def summarize(objects, centre, images):
    for obj in objects:
        pos = (obj['ra'],obj['dec'])
        obj['centre_distance_arcsec'] = separation(pos,centre)
        obj['image_distances_arcsec'] = [separation(pos,image) for image in images]
    objects = [obj for obj in objects if obj['centre_distance_arcsec'] <= CONE]
    # PS1 overlaps can return multiple rows for one objID. Same object can
    # never count as two distinct image counterparts. Resolve with catalog flags.
    objects.sort(key=lambda obj: (obj.get('primary_detection')!=1,obj.get('best_detection')!=1,
                                  -len(obj.get('valid_bands',[])),obj['centre_distance_arcsec']))
    unique = {}
    for obj in objects:
        unique.setdefault(obj['object_id'],obj)
    objects = sorted(unique.values(),key=lambda obj:obj['centre_distance_arcsec'])
    extended = [obj for obj in objects if obj['looks_extended']]
    if not extended:
        return dict(status='no_catalogued_extended_object' if objects else 'empty_catalogue_cone',
                    nearest_extended=None,objects=objects)
    near=extended[0]
    others=[other for other in objects if other['object_id']!=near['object_id']]
    assigned=[]
    if len(others)>=len(images)>0:
        distances=np.array([item['image_distances_arcsec'] for item in others]).T
        # Prefer assignments with the largest number of images inside 0.4".
        rows,indices=linear_sum_assignment(distances + (distances>IMAGE_RADIUS)*100)
        assigned=[dict(image_index=int(r),object_id=others[j]['object_id'],
                       distance_arcsec=float(distances[r,j])) for r,j in zip(rows,indices)]
    separated=(len(assigned)==len(images) and len(images)>=2 and
               all(x['distance_arcsec']<=IMAGE_RADIUS for x in assigned))
    delta=near['centre_distance_arcsec']
    closest=min(near['image_distances_arcsec'])
    central=delta<=CENTRE_RADIUS and delta+0.1<closest and closest>IMAGE_RADIUS
    if not central:
        status='extended_object_offset_or_on_image'
    elif not near['morphology_quality']:
        status='central_morphology_low_significance'
    elif not separated:
        status='central_extended_object_without_distinct_catalogued_images'
    else:
        status='central_extended_object_and_distinct_catalogued_images'
    return dict(status=status,nearest_extended=near,distinct_image_matches=assigned,
                central=central,all_images_separately_catalogued=separated,objects=objects)


def assess(systems,outdir):
    """Assess systems containing id, components with ra/dec, and center_radec_deg."""
    out=Path(outdir)
    out.mkdir(parents=True,exist_ok=True)
    summary=[]
    for system in systems:
        ident=str(system['id'])
        centre=list(map(float,system['center_radec_deg']))
        images=[(float(comp['ra']),float(comp['dec'])) for comp in system['components']]
        folder=out/ident
        folder.mkdir(exist_ok=True)
        evidence={}
        for survey in ('Legacy_DR10','Legacy_DR9','PS1_DR2_stack'):
            receipt={}
            try:
                rows,receipt=get_cone(survey,centre,folder)
                if rows is None:
                    result=dict(status='unavailable')
                else:
                    parser=panstarrs if survey=='PS1_DR2_stack' else tractor
                    objects=[obj for row in rows if (obj:=parser(row)) is not None]
                    result=summarize(objects,centre,images)
                    if not result.get('objects') and rows:
                        result['status']='no_usable_catalogue_rows'
            except Exception as exc:
                result=dict(status='unavailable',error=f'{type(exc).__name__}: {str(exc)[:400]}')
            evidence[survey]=dict(result=result,query=receipt)
            near=result.get('nearest_extended') or {}
            summary.append(dict(system_id=ident,survey=survey,status=result['status'],
                object_id=near.get('object_id'),ra=near.get('ra'),dec=near.get('dec'),
                centre_distance_arcsec=near.get('centre_distance_arcsec'),
                nearest_image_distance_arcsec=min(near['image_distances_arcsec']) if near else None))
        (folder/'catalogue_deblends.json').write_text(json.dumps(dict(centre_radec_deg=centre,
            image_source_ids=[str(row.get('source_id','')) for row in system['components']],
            surveys=evidence),indent=2,allow_nan=False)+'\n')
    labels=['system_id','survey','status','object_id','ra','dec',
            'centre_distance_arcsec','nearest_image_distance_arcsec']
    with (out/'catalogue_deblends.csv').open('w',newline='') as target:
        writer=csv.DictWriter(target,fieldnames=labels)
        writer.writeheader()
        writer.writerows(summary)
    return summary
