"""Matched-aperture Gaia FPR search. Execute on a GitHub Actions runner.

CDS I/361 provides both component positions and the accurate DR3 coordinates
of each parent. Duplicate observations in different parent fields are kept
through selection so the control and benchmark denominators remain explicit.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import NamedTuple
import gzip
import hashlib
import io
import json
import os

import numpy as np
import pandas as pd
import requests
from astropy.table import Table
from scipy.spatial import cKDTree

class Geometry(NamedTuple):
    name: str
    radius: float
    offset: float
    other_radius: float | None
    id_suffix: str

# Each field's assigned random bearing is identical between modes. The wide
# science and random footprints rotate into each other about their midpoint.
GEOMETRIES = {
    'primary': Geometry('primary',2.6,3.35,None,''),
    'wide': Geometry('wide',5.2,5.7,5.98,'W'),
}
MAX_IMAGE_SEPARATION = 5.2
SEED = 20231010
CHI2_LIMIT = 14.156   # 99.73% for three degrees of freedom
NATIVE = {
    'source_id': 'source_id', 'ra': 'ra', 'dec': 'dec', 'g': 'phot_g_mean_mag',
    'bp_rp': 'bp_rp', 'bp_snr': 'phot_bp_mean_flux_over_error',
    'rp_snr': 'phot_rp_mean_flux_over_error', 'ruwe': 'ruwe',
    'parallax': 'parallax', 'pmra': 'pmra', 'pmdec': 'pmdec',
    'parallax_error': 'parallax_error', 'pmra_error': 'pmra_error',
    'pmdec_error': 'pmdec_error', 'parallax_pmra_corr': 'parallax_pmra_corr',
    'parallax_pmdec_corr': 'parallax_pmdec_corr', 'pmra_pmdec_corr': 'pmra_pmdec_corr',
}
VIZIER = dict(zip(NATIVE, ['Source', 'RA_ICRS', 'DE_ICRS', 'Gmag', 'BP-RP',
    'RFBP', 'RFRP', 'RUWE', 'Plx', 'pmRA', 'pmDE', 'e_Plx', 'e_pmRA',
    'e_pmDE', 'PlxpmRAcor', 'PlxpmDEcor', 'pmRApmDEcor'], strict=True))
OUTCOLS = ['system_id','cohort','image_id','source_id','ra','dec','g','bp_rp',
           'field_id','input_lens_catalogue','fpr_ra','fpr_dec','fpr_g','n_obs','match_arcsec',
           'zero_chi2', *list(NATIVE)[5:]]
CHUNK = 2500
NATIVE_TAPS = [
    ('https://gea.esac.esa.int/tap-server/tap/sync','gaiadr3.gaia_source'),
    ('https://gaia.aip.de/tap/sync','gaiadr3.gaia_source'),
    ('https://irsa.ipac.caltech.edu/TAP/sync','gaia_dr3_source'),
    ('https://datalab.noirlab.edu/tap/sync','gaia_dr3.gaia_source'),
]


def unit_vectors(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])


def separation(ra1, dec1, ra2, dec2):
    """Exact spherical separation in arcseconds."""
    r1, d1, r2, d2 = np.radians([ra1, dec1, ra2, dec2])
    h = np.sin((d1-d2)/2)**2 + np.cos(d1)*np.cos(d2)*np.sin((r1-r2)/2)**2
    return np.degrees(2*np.arcsin(np.sqrt(np.clip(h, 0, 1))))*3600


def offset_points(ra, dec, offset):
    rng = np.random.default_rng(SEED)
    pa = rng.uniform(0, 2*np.pi, len(ra))
    dist = np.radians(offset/3600)
    ra, dec = np.radians(ra), np.radians(dec)
    d2 = np.arcsin(np.sin(dec)*np.cos(dist) + np.cos(dec)*np.sin(dist)*np.cos(pa))
    r2 = ra + np.arctan2(np.sin(pa)*np.sin(dist)*np.cos(dec),
                          np.cos(dist)-np.sin(dec)*np.sin(d2))
    return np.degrees(r2)%360, np.degrees(d2)


def read_fpr(cache, geometry):
    """Read fixed byte ranges from the I/361 dictionary, checking both row counts."""
    def parents():
        with gzip.open(cache/'lenscatn.dat.gz','rb') as stream:
            for line in stream:
                yield (int(line[:19]), int(line[40:42]), float(line[43:66]), float(line[67:90]))
    cats = np.fromiter(parents(), dtype=[('field_id','i8'),('cat','i1'),('ra','f8'),('dec','f8')])
    if len(cats) != 3760480:
        raise ValueError(f"Unexpected I/361/lenscatn row count {len(cats)}")
    counts = Counter()
    ngp = unit_vectors([192.85948], [27.12825])[0]

    def components():
        with gzip.open(cache/'lenscand.dat.gz','rb') as stream:
            for line in stream:
                n = int(line[71:74])
                counts['cds_component_rows'] += 1
                comp = int(line[75:78])
                if comp == 1:
                    counts['cds_fields_with_components'] += 1
                if not 2 <= n <= 8:
                    continue
                yield (int(line[:19]), comp, n, int(line[79:82]), int(line[83:85]),
                       float(line[86:109]), float(line[123:146]), float(line[202:212]),
                       float(line[226:249]), float(line[250:273]))
    obs = pd.DataFrame(np.fromiter(components(), dtype=[('field_id','i8'),('image_id','i2'),
        ('n_parent','i2'),('n_obs','i2'),('flag','i1'),('fpr_ra','f8'),('fpr_dec','f8'),
        ('fpr_g','f8'),('ra','f8'),('dec','f8')]))
    if counts['cds_component_rows'] != 4760920:
        raise ValueError(f"Unexpected I/361/lenscand row count: {counts}")
    parents = obs.loc[obs.image_id.eq(1), ['field_id','ra','dec','n_parent']].sort_values('field_id')
    if len(parents) != obs.field_id.nunique() or len(obs) != parents.n_parent.sum():
        raise ValueError("Incomplete multiplet groups in I/361")
    counts['fields_2_to_8_components'] = len(parents)
    vec = unit_vectors(parents.ra, parents.dec)
    parents = parents.loc[abs(vec@ngp) >= np.sin(np.radians(20))].copy()
    parents['input_lens_catalogue'] = parents.field_id.isin(cats['field_id'][np.isin(cats['cat'], [2,3])])
    obs = obs[obs.field_id.isin(parents.field_id)].drop(columns=['ra','dec'])
    obs['rid'] = np.arange(len(obs),dtype='int64')
    counts['science_fields'] = len(parents)
    control = parents.copy()
    control['ra'], control['dec'] = offset_points(parents.ra, parents.dec, geometry.offset)
    # A seed outside the intersection cannot contaminate a wide control. Check
    # all seeds in its first disk, because the nearest can lie outside disk two.
    tree = cKDTree(unit_vectors(cats['ra'],cats['dec']))
    control_xyz = unit_vectors(control.ra, control.dec)
    max_chord = 2*np.sin(np.radians(geometry.radius/3600)/2)
    if geometry.other_radius is None:
        good = tree.query(control_xyz, k=1)[0] > max_chord
    else:
        control['other_ra'], control['other_dec'] = parents.ra, parents.dec
        parents['other_ra'], parents['other_dec'] = control.ra, control.dec
        nearby = tree.query_ball_point(control_xyz, max_chord)
        sizes = np.fromiter((len(group) for group in nearby),dtype=int,count=len(nearby))
        owners = np.repeat(np.arange(len(nearby)),sizes)
        indices = np.fromiter((j for group in nearby for j in group),dtype=int,count=len(owners))
        inside = separation(parents.ra.to_numpy()[owners], parents.dec.to_numpy()[owners],
                            cats['ra'][indices], cats['dec'][indices]) <= geometry.other_radius
        good = np.ones(len(control),dtype=bool)
        good[owners[inside]]=False
    counts['control_apertures_on_fpr_qso'] = int((~good).sum())
    control = control.loc[good].copy()
    counts['control_fields'] = len(control)
    parents['cohort'], control['cohort'] = 'science','control'
    parents['system_id'] = 'F'+geometry.id_suffix+parents.field_id.astype(str)
    control['system_id'] = 'R'+geometry.id_suffix+control.field_id.astype(str)
    fields = pd.concat([parents,control],ignore_index=True)
    fields['stage'] = 'fewer_than_two_clean_fpr_components_in_aperture'
    return fields, obs, counts


def make_pairs(fields, obs, geometry):
    clean = obs[(obs.n_obs>=5) & obs.flag.isin([0,10]) & (obs.fpr_g<=21.0)]
    columns = ['system_id','field_id','ra','dec']
    if geometry.other_radius is not None:
        columns += ['other_ra','other_dec']
    members = fields[columns].merge(clean,on='field_id')
    r = separation(members.ra,members.dec,members.fpr_ra,members.fpr_dec)
    inside = r <= geometry.radius
    if geometry.other_radius is not None:
        r_other = separation(members.other_ra,members.other_dec,members.fpr_ra,members.fpr_dec)
        inside &= r_other <= geometry.other_radius
    members = members.loc[inside, ['system_id','rid','fpr_ra','fpr_dec','fpr_g']]
    p = members.merge(members,on='system_id',suffixes=('_a','_b'))
    p = p.loc[p.rid_a<p.rid_b].copy()
    fields.loc[fields.system_id.isin(p.system_id),'stage'] = 'pair_geometry_or_fpr_flux_contrast'
    p['sep_arcsec'] = separation(p.fpr_ra_a,p.fpr_dec_a,p.fpr_ra_b,p.fpr_dec_b)
    p = p.loc[p.sep_arcsec.between(0.70, MAX_IMAGE_SEPARATION) & ((p.fpr_g_a-p.fpr_g_b).abs()<=3.0)]
    p = p[['system_id','rid_a','rid_b','sep_arcsec']].copy()
    p['decision'] = 'no_unique_dr3_match'
    return p


def cds_match(sub, cache, receipts):
    """All positional associations first; recheck against DR3 epoch-2016 coordinates."""
    url='https://cdsxmatch.cds.unistra.fr/xmatch/api/v1/sync'
    payload=sub[['rid','fpr_ra','fpr_dec']].to_csv(index=False).encode()
    form={'request':'xmatch','distMaxArcsec':'2','RESPONSEFORMAT':'csv',
          'cat2':'vizier:I/355/gaiadr3','colRA1':'fpr_ra','colDec1':'fpr_dec','area':'allsky',
          'selection':'all'}
    key=hashlib.sha256(payload+json.dumps(form,sort_keys=True).encode()).hexdigest()
    dest, meta=cache/('xmatch-'+key+'.csv.gz'), cache/('xmatch-'+key+'.json')
    if dest.exists() and meta.exists():
        record=json.loads(meta.read_text())
        if record['cached_csv_gz_sha256']==hashlib.sha256(dest.read_bytes()).hexdigest():
            receipts.append(record)
            return pd.read_csv(dest,dtype={'source_id':str})
    response=requests.post(url,data=form,files={'cat1':('targets.csv',payload,'text/csv')},
                           timeout=(12,65))
    response.raise_for_status()
    data=pd.read_csv(io.BytesIO(response.content),dtype=str)
    names={v.lower():v for v in data}
    if 'rid' not in names or 'source' not in names:
        raise ValueError(f'CDS XMatch catalogue keys absent: {list(data)[:20]}')
    # The service may publish only a subset of the Gaia metadata. Fetch the
    # full covariance and photometry by exact (string-preserved) source IDs.
    keys=data[[names['rid'],names['source']]].copy()
    keys.columns=['rid','source_id']
    keys['rid']=pd.to_numeric(keys.rid,errors='raise').astype('int64')
    if len(keys)>=100000 or not keys.rid.isin(sub.rid).all() or not keys.source_id.str.fullmatch(r'[0-9]+').all():
        raise ValueError('Invalid CDS XMatch count, row ID, or source ID')
    native={k:names.get(v.lower()) for k,v in VIZIER.items()}
    native['ra']=native['ra'] or names.get('radeg')
    native['dec']=native['dec'] or names.get('dedeg')
    id_queries=[]
    if all(native.values()):
        df=keys.copy()
        for k,v in native.items():
            if k!='source_id':
                df[k]=pd.to_numeric(data[v],errors='raise')
    elif len(keys):
        blocks=[]
        ids=keys.source_id.drop_duplicates().tolist()
        for i in range(0,len(ids),700):
            block=ids[i:i+700]
            cols=', '.join(f'g.{v} AS {k}' for k,v in NATIVE.items())
            errors=[]
            for tap, table in NATIVE_TAPS:
                query=("SELECT "+cols+f" FROM {table} AS g WHERE g.source_id IN ("+
                       ','.join(block)+")")
                try:
                    resp=requests.post(tap,data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv',
                        'MAXREC':'2000','QUERY':query},timeout=(12,45))
                    resp.raise_for_status()
                    prop=pd.read_csv(io.BytesIO(resp.content),dtype=str)
                    prop.columns=prop.columns.str.lower()
                    if set(prop.columns)!=set(NATIVE) or set(prop.source_id)!=set(block) or len(prop)!=len(block):
                        raise ValueError('Incomplete source ID enrichment')
                    for k in NATIVE:
                        if k!='source_id':
                            prop[k]=pd.to_numeric(prop[k],errors='raise')
                    blocks.append(prop)
                    id_queries.append(dict(url=tap, query_sha256=hashlib.sha256(query.encode()).hexdigest(),
                        response_sha256=hashlib.sha256(resp.content).hexdigest(),n_ids=len(block)))
                    break
                except Exception as exc:
                    errors.append(str(exc)[:350])
            else:
                raise ValueError(f'CDS XMatch ID enrichment failed: {errors}')
        df=keys.merge(pd.concat(blocks,ignore_index=True),on='source_id',validate='many_to_one')
    else:
        df=pd.DataFrame(columns=['rid',*NATIVE])
    df.to_csv(dest,index=False)
    record=dict(url=url,parameters=form,n_upload=len(sub),n_returned=len(df),
                upload_sha256=hashlib.sha256(payload).hexdigest(),
                response_sha256=hashlib.sha256(response.content).hexdigest(),
                cached_csv_gz_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                id_enrichment=id_queries,first_rid=int(sub.rid.min()),last_rid=int(sub.rid.max()))
    meta.write_text(json.dumps(record))
    receipts.append(record)
    return df


def query_dr3(sub, cache, out, receipts, preferred):
    payload = io.BytesIO()
    Table.from_pandas(sub[['rid','fpr_ra','fpr_dec']]).write(payload, format='votable')
    upload = payload.getvalue()
    endpoints = [(url, table, NATIVE) for url, table in NATIVE_TAPS]
    endpoints.append(('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync','"I/355/gaiadr3"',VIZIER))
    endpoints.sort(key=lambda ep: ep[0]!=preferred[0])
    errors = []
    def try_cds():
        try:
            frame = cds_match(sub, cache, receipts)
            preferred[0] = 'cds'
            return frame
        except Exception as exc:
            errors.append(dict(url='CDS XMatch / DR3 by ID', error=str(exc)[:700]))
            print('CDS XMatch failed', errors[-1], flush=True)
            return None

    if preferred[0] in (None, 'cds'):
        frame = try_cds()
        if frame is not None:
            return frame
    for url, table, names in endpoints:
        def col(k):
            return 'g."'+names[k]+'"' if names is VIZIER else 'g.'+names[k]
        columns = ', '.join(f'{col(k)} AS {k}' for k in NATIVE)
        query = (f"SELECT p.rid AS rid, {columns} FROM TAP_UPLOAD.targets AS p JOIN {table} AS g "
                 f"ON 1=CONTAINS(POINT('ICRS',{col('ra')},{col('dec')}),"
                 "CIRCLE('ICRS',p.fpr_ra,p.fpr_dec,0.0000833333333333333))")
        if table == 'gaia_dr3_source':
            # IRSA TAP documents the cross-product form for uploaded table aliases.
            query = (f"SELECT p.rid AS rid, {columns} FROM TAP_UPLOAD.targets AS p, {table} AS g "
                     "WHERE CONTAINS(POINT(g.ra,g.dec),CIRCLE(p.fpr_ra,p.fpr_dec,"
                     "0.0000833333333333333))=1")
        key = hashlib.sha256(upload+query.encode()+url.encode()).hexdigest()
        dest = cache/('dr3-'+key+'.csv.gz')
        meta = cache/('dr3-'+key+'.json')
        try:
            if dest.exists() and meta.exists():
                record = json.loads(meta.read_text())
                if hashlib.sha256(dest.read_bytes()).hexdigest() != record['cached_csv_gz_sha256']:
                    raise ValueError("Cached extraction checksum mismatch")
                df = pd.read_csv(dest, dtype={'source_id':str})
            else:
                response = requests.post(url,
                    data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'votable','MAXREC':'100000',
                          'UPLOAD':'targets,param:targets','QUERY':query},
                    files={'targets':('targets.xml',upload,'application/x-votable+xml')},
                    timeout=(12,70))
                response.raise_for_status()
                import xml.etree.ElementTree as ET
                root = ET.fromstring(response.content)
                info = [(e.get('value','').upper(), e.text) for e in root.iter()
                        if e.tag.endswith('INFO') and e.get('name','').upper()=='QUERY_STATUS']
                if not info or any(value!='OK' for value,_ in info):
                    raise ValueError(f'TAP status {info!r}')
                df = Table.read(io.BytesIO(response.content), format='votable',
                                use_names_over_ids=True).to_pandas()
                df.columns = [s.lower() for s in df.columns]
                if set(df.columns) != {'rid', *NATIVE}:
                    raise ValueError(f'Incomplete DR3 schema {df.columns.tolist()}')
                if len(df) >= 100000 or not df.rid.isin(sub.rid).all():
                    raise ValueError('Invalid DR3 match count or row IDs')
                df['source_id'] = df.source_id.astype(str)
                if not df.source_id.str.fullmatch(r'[0-9]+').all():
                    raise ValueError('DR3 identifier lost integer precision')
                df.to_csv(dest, index=False)
                record = dict(url=url,query=query,upload_sha256=hashlib.sha256(upload).hexdigest(),
                              response_sha256=hashlib.sha256(response.content).hexdigest(),
                              cached_csv_gz_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                              n_upload=len(sub),n_returned=len(df),errors=errors)
                meta.write_text(json.dumps(record))
            receipts.append(record)
            preferred[0]=url
            return df
        except Exception as exc:
            errors.append(dict(url=url, error=f'{type(exc).__name__}: {str(exc)[:700]}'))
            print('DR3 provider failed', errors[-1], flush=True)
    if preferred[0] not in (None, 'cds'):
        frame = try_cds()
        if frame is not None:
            return frame
    raise RuntimeError(f'No complete DR3 response for {len(sub)} inputs: {errors}')


def covariance(row):
    # A 0.03 mas parallax and 0.03 mas/yr motion floor retain modest systematics.
    err = np.array([row.parallax_error,row.pmra_error,row.pmdec_error])
    corr = np.array([[1,row.parallax_pmra_corr,row.parallax_pmdec_corr],
                     [row.parallax_pmra_corr,1,row.pmra_pmdec_corr],
                     [row.parallax_pmdec_corr,row.pmra_pmdec_corr,1]])
    if not np.isfinite(corr).all() or np.any(abs(corr)>1) or not np.isfinite(err).all() or min(err)<=0:
        return None
    cov = np.outer(err,err)*corr + np.eye(3)*0.03**2
    return cov if np.linalg.eigvalsh(cov).min()>0 else None


def astrometric_selection(dr3):
    stats = {}
    for row in dr3.itertuples(index=False):
        c = covariance(row)
        # Broad DR3 ensemble zero point; the floor covers spatial/calibration variation.
        v = np.array([row.parallax+0.017,row.pmra,row.pmdec])
        quality = (c is not None and np.isfinite(v).all() and
                   0<row.ruwe<=1.6 and row.parallax_error<=1.5 and
                   max(row.pmra_error,row.pmdec_error)<=2.0)
        chi2 = float(v@np.linalg.solve(c,v)) if quality else np.nan
        phot = (np.isfinite(row.bp_rp) and min(row.bp_snr,row.rp_snr)>=5 and
                np.isfinite([row.bp_snr,row.rp_snr,row.g]).all())
        stats[row.rid] = (v,c,chi2,phot,row)
    return stats


def run(outdir):
    cache=Path(os.environ['GAIA_CACHE'])
    out=Path(outdir)
    out.mkdir(parents=True,exist_ok=True)
    mode = os.environ.get('GAIA_GEOMETRY','primary')
    if mode not in GEOMETRIES:
        raise ValueError(f'Choose GAIA_GEOMETRY from {tuple(GEOMETRIES)}')
    geometry = GEOMETRIES[mode]
    fields, obs, counts=read_fpr(cache,geometry)
    pairs=make_pairs(fields,obs,geometry)
    metadata = dict(name=geometry.name, radius_arcsec=geometry.radius,
                    other_radius_arcsec=geometry.other_radius, offset_arcsec=geometry.offset,
                    seed=SEED, max_image_separation_arcsec=MAX_IMAGE_SEPARATION)

    wanted=set(pairs.rid_a)|set(pairs.rid_b)
    targets=obs[obs.rid.isin(wanted)].sort_values('rid')
    print('FPR counts:',dict(counts),'DR3 positions:',len(targets),flush=True)
    fields.loc[fields.system_id.isin(pairs.system_id),'stage']='awaiting_dr3'
    fields.to_csv(out/'gaia_fields.csv.gz',index=False,float_format='%.9f')
    counts['dr3_target_count']=len(targets)
    (out/'gaia-fpr-selection.json').write_text(json.dumps(dict(counts=counts,geometry=metadata,
        fpr_pair_counts=pairs.groupby(pairs.system_id.str[0]).size().to_dict(),
        field_outcomes=fields.groupby(['cohort','stage']).size().unstack(fill_value=0).to_dict(orient='index')),
        indent=2)+'\n')
    matched, receipts, preferred, unavailable=[], [], [None], set()
    consecutive_failures=0
    for start in range(0,len(targets),CHUNK):
        batch=targets.iloc[start:start+CHUNK]
        if consecutive_failures>=2:
            unavailable.update(targets.rid.iloc[start:])
            receipts.append(dict(status='service_outage_deferred',n_inputs=len(targets)-start,
                                 first_rid=int(batch.rid.min())))
            break
        try:
            matched.append(query_dr3(batch,cache,out,receipts,preferred))
            consecutive_failures=0
        except Exception as exc:
            unavailable.update(batch.rid)
            receipts.append(dict(status='failed',first_rid=int(batch.rid.min()),
                last_rid=int(batch.rid.max()),n_inputs=len(batch),error=str(exc)[:3000]))
            consecutive_failures+=1
        (out/'gaia-dr3-queries.json').write_text(json.dumps(receipts,indent=2)+'\n')
        print('DR3 processed', min(len(targets),start+CHUNK),'/',len(targets),
              'unavailable:',len(unavailable), flush=True)
    (out/'gaia-dr3-queries.json').write_text(json.dumps(receipts,indent=2)+'\n')
    raw=pd.concat(matched,ignore_index=True) if matched else pd.DataFrame(columns=['rid',*NATIVE])
    raw=raw.merge(targets[['rid','fpr_ra','fpr_dec']],on='rid')
    raw['match_arcsec']=separation(*(raw[k].to_numpy(dtype=float) for k in ['ra','dec','fpr_ra','fpr_dec']))
    raw=raw[raw.match_arcsec<=0.30]
    # Multiple possible DR3 objects inside the FPR error circle are ambiguous.
    unique=raw.drop_duplicates('rid',keep=False)
    counts['dr3_component_queries_complete']=len(targets)-len(unavailable)
    counts['dr3_component_queries_unavailable']=len(unavailable)
    counts['unique_dr3_components']=len(unique)
    stats=astrometric_selection(unique)
    rows = []
    levels=['no_unique_dr3_match','same_dr3_source','missing_or_low_quality_astrometry',
            'inconsistent_with_zero_astrometry','differential_astrometry','missing_or_low_quality_colour',
            'colour_mismatch','passed_gaia']
    rank={k:i for i,k in enumerate(levels)}
    pair_results, peak, pairs_pass=[], {}, {}
    for p in pairs.itertuples(index=False):
        decision=levels[0]
        dchi=dc=np.nan
        a,b=stats.get(p.rid_a),stats.get(p.rid_b)
        if a is not None and b is not None:
            if a[4].source_id==b[4].source_id:
                decision=levels[1]
            elif not np.isfinite([a[2],b[2]]).all():
                decision=levels[2]
            elif max(a[2],b[2])>CHI2_LIMIT:
                decision=levels[3]
            else:
                v=a[0]-b[0]
                dchi=float(v@np.linalg.solve(a[1]+b[1],v))
                if dchi>CHI2_LIMIT:
                    decision=levels[4]
                elif not (a[3] and b[3]):
                    decision=levels[5]
                else:
                    dc=abs(a[4].bp_rp-b[4].bp_rp)
                    decision=levels[6] if dc>0.35 else levels[7]
        if decision=='passed_gaia':
            pairs_pass.setdefault(p.system_id,[]).append((p.rid_a,p.rid_b))
        if p.rid_a in unavailable or p.rid_b in unavailable:
            decision='dr3_service_unavailable'
        peak[p.system_id]=max(rank.get(decision,-1),peak.get(p.system_id,-1))
        pair_results.append((decision,dchi,dc))
    if pair_results:
        pairs[['decision','differential_chi2','abs_delta_bp_rp']]=pd.DataFrame(pair_results,index=pairs.index)
    else:
        pairs['differential_chi2']=np.nan
        pairs['abs_delta_bp_rp']=np.nan
    fields['stage']=fields.system_id.map({k:levels[v] if v>=0 else 'dr3_service_unavailable' for k,v in peak.items()}).fillna(fields.stage)
    uncertain=pairs.loc[pairs.decision.eq('dr3_service_unavailable'),'system_id']
    fields.loc[fields.system_id.isin(uncertain)&fields.stage.ne('passed_gaia'),'stage']='dr3_service_unavailable'
    fields['dr3_incomplete']=fields.system_id.isin(uncertain)
    stageindex=fields.set_index('system_id')
    obsindex=obs.set_index('rid')
    for system, edges in pairs_pass.items():
        # Use the largest mutually compatible image subset. Resolve ties by the
        # smallest total absolute colour difference, then FPR row order.
        nodes=sorted({x for edge in edges for x in edge})
        edges=set(map(lambda pair:tuple(sorted(pair)),edges))
        for size in range(len(nodes),1,-1):
            viable=[subset for subset in combinations(nodes,size)
                    if all(pair in edges for pair in combinations(subset,2))]
            if viable:
                best=min(viable,key=lambda subset:(sum(abs(stats[i][4].bp_rp-stats[j][4].bp_rp)
                                                    for i,j in combinations(subset,2)),subset))
                break
        ref=stageindex.loc[system]
        for i in best:
            native=stats[i][4]
            f=obsindex.loc[i]
            obj={k:getattr(native,k) for k in NATIVE}
            obj.update(system_id=system, cohort=ref.cohort, image_id=int(f.image_id),
                       field_id=str(ref.field_id),input_lens_catalogue=bool(ref.input_lens_catalogue),
                       fpr_ra=f.fpr_ra,fpr_dec=f.fpr_dec,fpr_g=f.fpr_g,n_obs=int(f.n_obs),
                       match_arcsec=native.match_arcsec,zero_chi2=stats[i][2])
            rows.append(obj)
    result=pd.DataFrame(rows,columns=OUTCOLS)
    result.to_csv(out/'gaia_components.csv',index=False)
    fields.to_csv(out/'gaia_fields.csv.gz', index=False, float_format='%.9f')
    pairs.to_csv(out/'gaia_pairs.csv.gz', index=False,float_format='%.6f')
    # Preserve all queried inputs, including missing measurements and failed cuts.
    audited=targets.merge(unique.drop(columns=['fpr_ra','fpr_dec']),how='left',on='rid')
    audited['dr3_coverage']=np.where(audited.rid.isin(unavailable),'unavailable','queried')
    audited['n_dr3_matches']=audited.rid.map(raw.groupby('rid').size()).fillna(0).astype('Int64')
    audited.loc[audited.rid.isin(unavailable),'n_dr3_matches']=pd.NA
    audited['zero_chi2']=audited.rid.map({k:v[2] for k,v in stats.items()})
    audited.to_csv(out/'gaia_dr3_audit.csv.gz',index=False,float_format='%.10g')
    totals=fields.groupby(['cohort','stage']).size().unstack(fill_value=0).to_dict(orient='index')
    (out/'gaia-counts.json').write_text(json.dumps(dict(counts=counts,geometry=metadata,field_outcomes=totals,
        controls='Random offset in each eligible FPR field; reject a control containing any FPR input in its full selection footprint. Identical FPR/DR3 cuts apply. The random-position rate is conditional on selected FPR environments.', 
        schema='CDS I/361, frozen 2023 files; not a count of the live AIP mirror.'),indent=2)+'\n')
    return result, fields


if __name__=='__main__':
    import sys
    run(sys.argv[1])
