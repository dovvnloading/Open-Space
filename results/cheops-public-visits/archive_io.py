"""Anonymous CHEOPS inventory and products. Execute network operations on the runner.

DACE routes follow astro-dace/dace-query, src/dace_query/cheops/cheops.py.
ESASky exposes the mission archive's product links in observations.mv_cheops_obs_fdw.
"""
import argparse
import configparser
import csv
import gzip
import hashlib
import io
import json
import math
import re
import tarfile
import tempfile
from datetime import datetime, timezone
from importlib import metadata
from functools import lru_cache
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

import requests
from astropy.io import fits

TAP = 'https://sky.esa.int/esasky-tap/tap/sync'
TABLE = 'observations.mv_cheops_obs_fdw'
PORTALS = ('https://cheops.unige.ch/archive_browser/',
           'https://cheops-archive.ssdc.asi.it/archive_browser/')
LIMIT = 100000
PIPE_NAME = re.compile(r'(?:^|[/_ -])(?:pipe|psf)(?:$|[/_. -])|_(?:sa|im)\.fits(?:\.gz)?$', re.I)


@lru_cache(maxsize=1)
def dace_host():
    try:
        p = metadata.distribution('dace-query').locate_file('dace_query/config.ini')
        cfg = configparser.ConfigParser()
        cfg.read(p)
        return cfg['api']['cheops-webapp'].rstrip('/')
    except (KeyError, metadata.PackageNotFoundError):
        return 'https://cheops-webapp-pg.obsuksprd2.unige.ch'


def log(path, **values):
    with Path(path).open('a', encoding='utf-8') as out:
        out.write(json.dumps(values, sort_keys=True, allow_nan=False) + '\n')


def digest(path):
    with Path(path).open('rb') as src:
        h = hashlib.sha256()
        for block in iter(lambda: src.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def exchange(evidence, label, url, *, params=None, body=None, local=None):
    """Hash the exact decoded response bytes before interpreting any successful payload."""
    row = dict(operation=label, url=url, utc=datetime.now(timezone.utc).isoformat())
    if params is not None:
        row['query'] = params
    if body is not None:
        row['body'] = body
    try:
        with requests.Session() as session:
            # No .netrc, ambient auth, or archive API keys.
            session.trust_env = False
            session.headers.update({'User-Agent': 'Open-Space CHEOPS archive research/1.0',
                                    'Accept': ('application/octet-stream' if local else
                                               'application/json' if label.startswith('DACE') else '*/*')})
            with session.request('POST' if body is not None else 'GET', url, params=params,
                                 json=body, timeout=(10, 100), stream=True) as response:
                row.update(http=response.status_code, content_type=response.headers.get('Content-Type'),
                           final_url=response.url.partition('?')[0])
                if response.status_code >= 400:
                    sample = next(response.iter_content(2048), b'')
                    row['response_excerpt'] = re.sub(r'\s+', ' ', sample.decode(errors='replace'))[:450]
                    row['auth_challenge'] = bool(response.status_code == 401 or
                        response.headers.get('WWW-Authenticate') or (response.status_code == 403 and
                        re.search(r'authentication required|please (?:log|sign) in|missing (?:api key|credentials)',
                                  row['response_excerpt'], re.I)))
                    response.raise_for_status()
                out = local.open('wb') if local else io.BytesIO()
                with out:
                    h = hashlib.sha256()
                    size = 0
                    for buf in response.iter_content(1024 * 1024):
                        h.update(buf)
                        size += len(buf)
                        out.write(buf)
                    row.update(sha256=h.hexdigest(), bytes=size)
                    if local is None:
                        return out.getvalue()
                    return {'sha256': h.hexdigest(), 'disposition':
                            response.headers.get('Content-Disposition', '')}
    except Exception as exc:
        row['error'] = f'{type(exc).__name__}: {str(exc)[:300]}'
        raise
    finally:
        log(evidence, **row)


def unpack_columns(data):
    """DACE serializes a response as protobuf columns with optional run lengths."""
    if 'parameters' not in data:
        raise ValueError('DACE response lacks parameter columns')
    cols = {}
    for p in data['parameters']:
        vals = next((p[k] for k in ('stringValues','intValues','longValues','doubleValues',
                                    'floatValues','boolValues') if k in p), [])
        def expand(values):
            if p.get('occurrences'):
                return [v for v,n in zip(values,p['occurrences']) for _ in range(n)]
            return values
        cols[p['variableName']] = expand(vals)
        if 'minErrorValues' in p:
            cols[p['variableName'] + '_err'] = expand(p['minErrorValues'])
    return cols


def number(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def canonical(key):
    m = re.search(r'(?:CH_)?(PR\d+_TG\d+)', str(key), re.I)
    return 'CH_' + m[1].upper() if m else None


def dace_inventory(evidence):
    cols = unpack_columns(json.loads(exchange(evidence, 'DACE visit listing', dace_host()+'/search',
                                   params={'limit':LIMIT,'sort':'{}','filters':'{}'})))
    keys = cols['file_key']
    if len(cols['status_published']) != len(keys):
        raise ValueError('missing per-visit public flag')
    best = {}
    for i,key in enumerate(keys):
        if cols['status_published'][i] is not True:
            continue
        def val(name):
            values = cols.get(name, [])
            return values[i] if i < len(values) else None
        stem = canonical(key) or str(key)
        if stem in best and str(best[stem]['file_key']) >= str(key):
            continue
        ra,dec = number(val('ra_deg')),number(val('dec_deg'))
        if (ra is None or dec is None) and val('obj_pos_coordinates_hms_dms'):
            from astropy.coordinates import SkyCoord
            from astropy import units as u
            try:
                pos = SkyCoord(val('obj_pos_coordinates_hms_dms'),unit=(u.hourangle,u.deg))
                ra,dec = float(pos.ra.deg),float(pos.dec.deg)
            except (TypeError,ValueError):
                pass
        best[stem] = dict(visit_id=stem, file_key=str(key), star=val('obj_id_catname'),
                          obs_id=str(val('obs_id')), ra=ra,
                          dec=dec, coords_hms_dms=val('obj_pos_coordinates_hms_dms'),
                          start_mjd=number(val('date_mjd_start')), pi=val('pi_name'),
                          file_rootpath=val('file_rootpath'),
                          dace=True, links=[])
    log(evidence, operation='DACE listing summary', rows=len(keys), visits=len(best),
        public=sum(v is True for v in cols['status_published']))
    return best, len(keys) < LIMIT


def tap(evidence, adql):
    raw = exchange(evidence, 'ESA TAP', TAP, params={'REQUEST':'doQuery', 'LANG':'ADQL',
                                'FORMAT':'csv','MAXREC':LIMIT,'QUERY':adql})
    if raw.lstrip().startswith(b'<'):
        raise ValueError('TAP returned XML or HTML: '+raw[-400:].decode(errors='replace'))
    return list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))


def esa_inventory(evidence, visits):
    rows = tap(evidence,'SELECT * FROM '+TABLE)
    total = tap(evidence,'SELECT COUNT(*) AS n FROM '+TABLE)
    ids = {r['obs_id']: k for k,r in visits.items() if r['obs_id'] not in ('','None')}
    for raw in rows:
        r = {k.lower():v for k,v in raw.items()}
        oid = r.get('obs_id','')
        links = []
        for k,v in r.items():
            if v and urlsplit(v).scheme in ('http','https') and re.search(
                    r'sci_cor_lc|lightcurve|_lc_|subarray|starcatalog|_pipe|psf|_(sa|im)\.fits',k+' '+v,re.I):
                links.append({'column':k,'url':v})
        stem = ids.get(oid) or canonical(' '.join(str(v) for v in r.values()))
        if not stem:
            if not oid:
                raise ValueError('ESA visit missing observation identifier')
            stem = 'ESA_'+oid
        if stem not in visits:
            visits[stem] = dict(visit_id=stem, file_key=stem, star=r.get('targname'), obs_id=oid,
                                ra=number(r.get('ra_deg')),dec=number(r.get('dec_deg')),
                                dace=False, links=[])
        visits[stem]['links'].extend(x for x in links if x not in visits[stem]['links'])
        for key in ('ra','dec'):
            if visits[stem].get(key) is None:
                visits[stem][key] = number(r.get(key+'_deg'))
    return len(rows) == int(next(iter(total[0].values())))


def safe_name(value):
    return re.sub(r'[^A-Za-z0-9_.+-]', '_', Path(unquote(value)).name)[:210] or 'product'


def keep_fits(payload, dest, suggested, evidence, transfer):
    """Accept FITS by signature, never by HTTP 200, and flatten tar members safely."""
    files = {}
    with payload.open('rb') as f:
        magic = f.read(8)
    compressed = magic[:2] == b'\x1f\x8b'
    if compressed:
        with gzip.open(payload,'rb') as f:
            magic = f.read(8)
    if magic == b'SIMPLE  ':
        name = re.search(r'filename="?([^";]+)',transfer['disposition'],re.I)
        if (name and name[1].lower().endswith(('.fits','.fits.gz')) and
                re.search(r'sci_|ext_|pipe|psf|_(sa|im)\.fits',name[1],re.I)):
            suggested = name[1]
        if not suggested.lower().endswith(('.fits','.fits.gz')):
            suggested += '.fits'
        path = dest/safe_name(suggested)
        if compressed and not path.name.endswith('.gz'):
            path = Path(str(path)+'.gz')
        payload.replace(path)
        files[path] = transfer['sha256']
    elif tarfile.is_tarfile(payload):
        with tarfile.open(payload,'r:*') as archive:
            for member in archive:
                if not member.isfile() or not member.name.lower().endswith(('.fits','.fits.gz')):
                    continue
                path = dest/safe_name(member.name)
                with archive.extractfile(member) as src, path.open('wb') as out:
                    h = hashlib.sha256()
                    for block in iter(lambda:src.read(1024 * 1024),b''):
                        h.update(block)
                        out.write(block)
                files[path] = h.hexdigest()
    else:
        with payload.open('rb') as f:
            sample = f.read(400).decode(errors='replace')
        raise ValueError('response was neither FITS nor an archive: '+sample[:160])
    valid = []
    for p,sha in files.items():
        with fits.open(p,memmap=False,lazy_load_hdus=True) as hdus:
            hdr = hdus[1].header if len(hdus)>1 else hdus[0].header
            log(evidence,operation='FITS',file=p.name,sha256=sha,bytes=p.stat().st_size,
                extname=hdr.get('EXTNAME'),object=hdr.get('OBJECT'),naxis2=hdr.get('NAXIS2'))
        valid.append(p)
    if not valid:
        raise ValueError('archive contained no FITS products')
    return valid


def download_product(record, dest, evidence, dtype, explicit_files=None, aperture=None):
    filters = {'file_key':{'equal':[record['file_key']]}}
    body = {'fileType':dtype,'filters':filters,'aperture':aperture}
    if explicit_files:
        body = {'fileType':dtype,'files':explicit_files}
    ticket = json.loads(exchange(evidence,'DACE request '+dtype,dace_host()+'/download',
                        body=body))['key']
    payload = dest/'download.tmp'
    transfer = exchange(evidence, 'DACE anonymous product '+dtype,
             dace_host()+'/download/photometry/'+quote(str(ticket),safe=''),
             params={'compressed':'false'},local=payload)
    postfix = '-'+aperture.upper() if aperture else ''
    title = record['file_key']+'_'+dtype+postfix+'.fits'
    if dtype=='files' and explicit_files and len(explicit_files)==1:
        title = Path(explicit_files[0]).name
    return keep_fits(payload,dest,title,evidence,transfer)


def dace_timeseries(record,dest,evidence):
    """Fallback for API users allowed vector photometry but denied packaged FITS."""
    import numpy as np
    filters = {'file_key':{'equal':[record['file_key']]}}
    path = dace_host()+'/photometry/'+quote(record['star'],safe='')
    args = dict(aperture='default', filters=json.dumps(filters),sort='{}')
    try:
        raw = exchange(evidence,'DACE DRP flux vectors',path,body=args)
    except requests.HTTPError as exc:
        if exc.response.status_code not in (404,405):
            raise
        raw = exchange(evidence,'DACE DRP flux vectors legacy',path,params=args)
    cols = unpack_columns(json.loads(raw))
    if not cols.get('obj_date_bjd_vect') or not cols.get('photom_flux_vect'):
        raise ValueError('no vector photometry for requested visit')
    if set(map(str,cols.get('file_key',[]))) - {record['file_key']}:
        raise ValueError('vector response contains multiple visits')
    names = {'obj_date_bjd_vect':'BJD_TIME','photom_flux_vect':'FLUX','photom_flux_vect_err':'FLUXERR',
             'photom_roll_angle_vect':'ROLL_ANGLE','photom_conta_lc_vect':'CONTA_LC',
             'photom_centroid_x_vect':'CENTROID_X','photom_centroid_y_vect':'CENTROID_Y',
             'photom_background_vect':'BACKGROUND','photom_smearing_lc_vect':'SMEARING_LC',
             'photom_event_vect':'EVENT'}
    columns = []
    for src,dst in names.items():
        if src not in cols:
            continue
        v = np.array(cols[src],dtype=float)
        if dst == 'EVENT':
            for name,vals in cols.items():
                if (name.startswith('photom_event_') and name.endswith('_flag')
                        and len(vals)==len(v)):
                    v[np.asarray(vals,dtype=bool)] = 1.
        if dst == 'BJD_TIME' and np.nanmedian(v)<100000:
            v += 2400000.
        columns.append(fits.Column(name=dst,format='D',array=v))
    hdu = fits.BinTableHDU.from_columns(columns)
    hdu.header['OBJECT'] = record['star']
    if record.get('ra') is not None and record.get('dec') is not None:
        hdu.header['RA_TARG'] = record['ra']
        hdu.header['DEC_TARG'] = record['dec']
    output = dest/(safe_name(record['file_key'])+'_DRP_DEFAULT_vector.fits')
    hdu.writeto(output,overwrite=True)
    log(evidence,operation='derived FITS from DACE vectors',file=output.name,
        sha256=digest(output),source_sha256=hashlib.sha256(raw).hexdigest())
    return [output]


def fetch_visit(record,dest,include_field=False):
    dest = Path(dest)
    dest.mkdir(parents=True,exist_ok=True)
    evidence = dest/'manifest.jsonl'
    files = []
    def attempt(operation,fn,*args):
        try:
            return fn(*args)
        except Exception as exc:
            log(evidence,operation=operation,error=f'{type(exc).__name__}: {str(exc)[:400]}')
            return []
        finally:
            (dest/'download.tmp').unlink(missing_ok=True)
    drp = []
    if record['dace']:
        for aperture in ('default','optimal'):
            drp = attempt('DRP '+aperture,download_product,record,dest,evidence,'lightcurves',None,aperture)
            if drp:
                break
        if not drp:
            drp = attempt('DRP aperture bundle',download_product,record,dest,evidence,'lightcurves')
        if not drp and record.get('file_rootpath'):
            drp = attempt('DRP legacy file-list route',download_product,record,dest,evidence,
                          'lightcurves',[record['file_rootpath']])
        files.extend(drp)
    for link in record['links']:
        name = link['column']+' '+link['url']
        is_field = bool(re.search('subarray|starcatalog',name,re.I))
        is_pipe = bool(PIPE_NAME.search(name))
        if (is_field and not include_field) or (not is_pipe and not is_field and drp):
            continue
        def get_link():
            payload = dest/'download.tmp'
            transfer = exchange(evidence,'ESA '+link['column'],link['url'],local=payload)
            filename = Path(urlsplit(link['url']).path).name
            if '.fits' not in filename.lower():
                filename = record['file_key']+'_'+link['column']+'.fits'
            return keep_fits(payload,dest,filename,evidence,transfer)
        found = attempt('ESA linked FITS',get_link)
        files.extend(found)
        if found and not is_field and not is_pipe:
            drp = found
    if not drp and record['dace']:
        files.extend(attempt('DACE vector fallback',dace_timeseries,record,dest,evidence))
    if include_field and record['dace']:
        listed = sample_products(evidence,record)
        for name in listed:
            if PIPE_NAME.search(name) and not any(f.name==Path(name).name for f in files):
                files.extend(attempt('catalogued PIPE file',download_product,record,dest,evidence,
                                     'files',[name]))
        catalogues = [p for p in files if 'starcatalog' in p.name.lower()]
        if not catalogues:
            catalogues = attempt('field catalogue',download_product,record,dest,evidence,'EXT_PRE_StarCatalogue')
            files.extend(catalogues)
        usable = False
        for cat in catalogues:
            try:
                with fits.open(cat,memmap=False) as hdus:
                    usable |= any({'RA','DEC'}.issubset(h.columns.names or []) and h.header.get('NAXIS2',0)>1
                                  for h in hdus if isinstance(h,fits.BinTableHDU))
            except (ValueError,OSError) as exc:
                log(evidence,operation='field catalogue validation',error=str(exc)[:350])
        if usable and not any('sci_cor_subarray' in f.name.lower() for f in files):
            files.extend(fetch_localization_images(record,dest))
        elif not usable:
            log(evidence,operation='corrected subarray',skip='no catalogue with multiple source coordinates')
    log(evidence,operation='visit product summary',visit_id=record['visit_id'],
        files=[f.name for f in files],
        pipe_identified=any(PIPE_NAME.search(f.name) is not None for f in files))
    return list(dict.fromkeys(files))


def fetch_localization_images(record,dest):
    """Get a correction-stage cube when a viable catalogue or a target event needs it."""
    dest = Path(dest)
    evidence = dest/'manifest.jsonl'
    existing = [p for p in dest.glob('*.fits*') if 'sci_cor_subarray' in p.name.lower()]
    if existing or not record['dace']:
        return existing
    try:
        return download_product(record,dest,evidence,'SCI_COR_SubArray')
    except Exception as exc:
        log(evidence,operation='corrected subarray',error=str(exc)[:350])
        return []
    finally:
        (dest/'download.tmp').unlink(missing_ok=True)


def fetch_validation_apertures(record,dest):
    """Independent aperture sizes for a detected target event; never count as separate visits."""
    dest = Path(dest)
    evidence = dest/'manifest.jsonl'
    result = []
    if not record['dace']:
        return result
    for aperture in ('rinf','rsup'):
        existing = sorted(p for p in dest.glob('*.fits*') if aperture in p.name.lower())
        if existing:
            result.extend(existing)
            continue
        try:
            result.extend(download_product(record,dest,evidence,'lightcurves',aperture=aperture))
        except Exception as exc:
            log(evidence,operation='aperture validation '+aperture,error=str(exc)[:350])
        finally:
            (dest/'download.tmp').unlink(missing_ok=True)
    return result


def sample_products(evidence,record):
    """Inspect the archive product catalogue for the verified observation."""
    field = re.compile(r'subarray|starcatalog', re.I)
    names = []
    summary = dict(operation='sample product inventory', file_key=record['file_key'],
                   ESA_columns=[x['column'] for x in record['links']])
    if record['dace']:
        body = {'fileType':'all','aperture':None,'filters':{'file_key':{'equal':[record['file_key']]}}}
        try:
            raw = exchange(evidence,'DACE browse all product names',dace_host()+'/download/browse',body=body)
            cols = unpack_columns(json.loads(raw))
            names = [str(v) for v in cols['file']]
            summary.update(listed_files=len(names), field_product_examples=[
                Path(n).name for n in names if field.search(n)][:10],
                PIPE_name_matches=[Path(n).name for n in names if PIPE_NAME.search(n)][:20])
        except Exception as exc:
            summary['listing_error'] = f'{type(exc).__name__}: {str(exc)[:350]}'
    log(evidence,**summary)
    return names


def inventory(out):
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    evidence = out/'archive_probes.jsonl'
    evidence.write_text('')
    rows = {}
    dace_complete = esa_complete = False
    try:
        rows,dace_complete = dace_inventory(evidence)
    except Exception as exc:
        log(evidence,operation='DACE inventory error',error=str(exc)[:500])
    try:
        esa_complete = esa_inventory(evidence,rows)
    except Exception as exc:
        log(evidence,operation='ESA inventory error',error=str(exc)[:500])
    if not rows:
        for url in PORTALS:
            try:
                exchange(evidence,'mission portal (landing page only)',url)
            except requests.RequestException:
                pass
    records = sorted(rows.values(),key=lambda r:r['visit_id'])
    with (out/'inventory.jsonl').open('w') as f:
        for r in records:
            f.write(json.dumps(r,sort_keys=True,allow_nan=False)+'\n')
    probes = [json.loads(line) for line in evidence.read_text().splitlines()]
    listing = [r for r in probes if r['operation'] in ('DACE visit listing','ESA TAP')]
    status = 'ready' if rows else 'unavailable'
    if not rows and len(listing)>=2 and all(r.get('auth_challenge') for r in listing):
        status='access_block'
    tested = []
    success = False
    probe_auth = []
    if records:
        # Choose visits distributed across programmes and delivery routes. A sample
        # can establish usability; its failure only bounds this access attempt.
        selections = []
        for group in ([r for r in records if r['dace'] and r['links']],
                      [r for r in records if r['dace']], [r for r in records if r['links']]):
            if group:
                for i in (len(group)//2, 0, len(group)-1):
                    if group[i] not in selections:
                        selections.append(group[i])
        for r in selections:
            tested.append(r['file_key'])
            with tempfile.TemporaryDirectory() as temp:
                products = fetch_visit(r, Path(temp))
                visit_log = (Path(temp)/'manifest.jsonl').read_text()
                with evidence.open('a') as outlog:
                    outlog.write(visit_log)
                attempt_log = [json.loads(line) for line in visit_log.splitlines()]
                direct = [v for v in attempt_log if v.get('url') and (
                    v['operation'].startswith(('DACE anonymous product','DACE DRP flux vectors','ESA ')) or
                    (v['operation'].startswith('DACE request') and v.get('error')))]
                probe_auth.append(bool(direct) and all(v.get('auth_challenge') for v in direct))
                for product in products:
                    with fits.open(product,memmap=False) as hdus:
                        for hdu in hdus:
                            names = {n.upper() for n in (getattr(getattr(hdu,'columns',None),'names',[]) or [])}
                            if (names.intersection({'FLUX','FLUX_PSF','PSF_FLUX'}) and
                                    names.intersection({'BJD_TIME','BJD','BJD_TDB','TIME'}) and
                                    hdu.header.get('NAXIS2',0) > 1):
                                import numpy as np
                                tcol = next(n for n in ('BJD_TIME','BJD','BJD_TDB','TIME') if n in names)
                                ts = np.asarray(hdu.data[tcol], dtype=float)
                                ts = ts[np.isfinite(ts)]
                                if not ts.size:
                                    continue
                                header_keys = ('TELESCOP','INSTRUME','EXTNAME','PROC_CHN','PIPE_VER',
                                               'TIMESYS','BJDREFI','BJDREFF','MJDREFI','TIMEUNIT',
                                               'TARGNAME','OBJECT','RA_TARG','DEC_TARG','AP_TYPE')
                                header = dict(hdus[0].header)
                                header.update(dict(hdu.header))
                                log(evidence,operation='verified anonymous light curve',
                                    file_key=r['file_key'],file=product.name,sha256=digest(product),
                                    rows=len(hdu.data),time_column=tcol,time_min=float(np.min(ts)),
                                    time_max=float(np.max(ts)),
                                    columns={c.name:str(c.unit or '') for c in hdu.columns},
                                    header={k:str(header[k]) for k in header_keys if k in header})
                                success = True
                if success:
                    sample_products(evidence,r)
                    break
        if not success:
            status = 'access_block' if probe_auth and all(probe_auth) else 'unavailable'
    result = dict(status=status, visits=len(records), tested_file_keys=tested,
                  verified_lightcurve=success,inventory_complete=bool(records and dace_complete),
                  dace_complete=dace_complete,esa_complete=esa_complete,
                  inventory_sha256=digest(out/'inventory.jsonl'))
    (out/'archive_status.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command',choices=['inventory'])
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    inventory(args.directory)


def fetch_repair_inputs(visit, old_record, dest):
    """Restore the three hash-pinned FITS products used in the original image attempt.

    The archived DACE download ticket can be reused while live. If it has
    expired, obtain a fresh ticket for the same visit, type and aperture.
    """
    dest = Path(dest)
    dest.mkdir(parents=True,exist_ok=True)
    evidence = dest/'manifest.jsonl'
    observation = old_record['visit_id']
    if observation != visit['visit_id']:
        raise ValueError('Visit and original receipts disagree')
    latest = None
    found = {}
    primary = next((p.get('name') for p in old_record.get('products',[])
                    if p.get('source_kind')=='target' and p.get('reduction')=='DRP'),None)
    if primary is None:
        primary = next((p.get('name') for p in old_record.get('products',[]) if p.get('reduction')=='DRP'),None)
    for row in old_record.get('receipts',[]):
        if row.get('operation','').startswith('DACE anonymous product'):
            latest = row
        if row.get('operation')!='FITS' or not row.get('sha256'):
            continue
        name=row['file']
        if name==primary:
            kind='lightcurves'
        elif 'starcatalog' in name.lower():
            kind='EXT_PRE_StarCatalogue'
        elif 'sci_cor_subarray' in name.lower():
            kind='SCI_COR_SubArray'
        else:
            continue
        if kind not in found:
            found[kind]=(name,row['sha256'],latest)
    types=('lightcurves','EXT_PRE_StarCatalogue','SCI_COR_SubArray')
    if set(found)!=set(types):
        raise ValueError('Missing original hash-pinned input(s): '+', '.join(k for k in types if k not in found))
    paths=[]
    for kind in types:
        name,expected,ticket=found[kind]
        payload=dest/'download.tmp'
        url=ticket.get('url','') if ticket else ''
        # Use precisely the published download URL and never attach auth.
        allowed=(urlsplit(url).scheme=='https' and urlsplit(url).hostname==
                 urlsplit(dace_host()).hostname and '/download/photometry/' in urlsplit(url).path)
        retrieved=[]
        receipt_start=evidence.stat().st_size if evidence.exists() else 0
        if allowed:
            try:
                transfer=exchange(evidence,'restore original '+kind,url,params=ticket.get('query'),local=payload)
                retrieved=keep_fits(payload,dest,name,evidence,transfer)
            except Exception as exc:
                log(evidence,operation='original ticket unavailable',file=name,error=str(exc)[:350])
            finally:
                payload.unlink(missing_ok=True)
        if not retrieved:
            aperture=None
            if kind=='lightcurves':
                match=re.search(r'-(DEFAULT|OPTIMAL|RINF|RSUP)[_.]',name,re.I)
                aperture=match[1].lower() if match else None
            try:
                retrieved=download_product(visit,dest,evidence,kind,aperture=aperture)
            finally:
                payload.unlink(missing_ok=True)
        # Reuse stream/extraction checksums rather than re-reading large cubes.
        with evidence.open('rb') as src:
            src.seek(receipt_start)
            verified={r['file'] for raw in src for r in (json.loads(raw),)
                      if r.get('operation')=='FITS' and r.get('sha256')==expected}
        match=next((p for p in retrieved if p.name in verified),None)
        if match is None:
            raise ValueError('Unable to reproduce original FITS checksum: '+name)
        log(evidence,operation='original input checksum verified',file=name,sha256=expected)
        paths.append(match)
    return paths, paths[0]
