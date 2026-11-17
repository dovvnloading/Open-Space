"""Download catalog neighbors and direct optical/WISE image cutouts at candidate positions."""
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
import sys
import time
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'wise'))
from fetch import targets


def get(url):
    for attempt in range(3):
        try:
            with urlopen(url, timeout=70) as response:
                return response.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def measure_catalogs(stars):
    clauses = [f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{float(t['ra']):.8f},{float(t['dec']):.8f},0.0041666667))"
               for t in stars]
    surveys = [
        ('gaia','https://gea.esac.esa.int/tap-server/tap/sync','gaiadr3.gaia_source',
         'source_id,ra,dec,pmra,pmdec,ruwe,phot_g_mean_mag,non_single_star'),
        ('simbad','https://simbad.cds.unistra.fr/simbad/sim-tap/sync','basic',
         'main_id,otype,ra,dec,plx_value')]
    receipts = []
    for name, endpoint, table, columns in surveys:
        query = f'SELECT {columns} FROM {table} WHERE ' + ' OR '.join(clauses)
        url = endpoint + '?' + urlencode(dict(LANG='ADQL', FORMAT='csv', REQUEST='doQuery', MAXREC=2000, QUERY=query))
        receipt = dict(source=name, utc=dt.datetime.now(dt.timezone.utc).isoformat(), url=url)
        try:
            body = get(url)
            rd = csv.DictReader(io.StringIO(body.decode('utf-8-sig')))
            if not {'ra','dec'}.issubset(rd.fieldnames or []):
                raise ValueError(body[:250])
            count = sum(1 for _ in rd)
            if count >= 2000:
                raise ValueError('MAXREC reached')
            (ROOT / (name + '_neighborhood.csv')).write_bytes(body)
            receipt.update(rows=count, sha256=hashlib.sha256(body).hexdigest())
        except Exception as error:
            receipt['error'] = str(error)
        receipts.append(receipt)
    (ROOT/'neighborhood_queries.json').write_text(json.dumps(receipts, indent=2)+'\n')


def cutout(task):
    star, layer = task
    args = dict(ra=star['ra'], dec=star['dec'], size=120, pixscale=.4, layer=layer)
    url = 'https://www.legacysurvey.org/viewer/jpeg-cutout?' + urlencode(args)
    receipt = dict(label=star['label'], layer=layer, url=url, utc=dt.datetime.now(dt.timezone.utc).isoformat())
    try:
        body = get(url)
        im = Image.open(io.BytesIO(body)).convert('RGB')
        if im.size != (120,120):
            raise ValueError(f'unexpected image size {im.size}')
        arr = np.asarray(im)
        status = 'displayed' if np.ptp(arr) else 'empty image'
        receipt.update(sha256=hashlib.sha256(body).hexdigest(), bytes=len(body), status=status)
        return receipt, im
    except Exception as error:
        receipt.update(status='unavailable', error=str(error)[:300])
        return receipt, Image.new('RGB',(120,120),'#111111')


def main():
    stars = targets()
    measure_catalogs(stars)
    tasks = [(t,layer) for t in stars for layer in ('ls-dr10','unwise-neo6')]
    with ThreadPoolExecutor(4) as pool:
        results = list(pool.map(cutout,tasks))
    sheet = Image.new('RGB',(544,762),'white')
    draw = ImageDraw.Draw(sheet)
    for i,(rec,im) in enumerate(results):
        j = i//2
        y, x = 21+146*(j%5), 10+132*(2*(j//5)+i%2)
        draw.text((x,y-15),rec['label']+'  '+rec['layer'],fill='#111111')
        sheet.paste(im,(x,y))
        if rec['status'] != 'displayed':
            draw.text((x+4,y+50),rec['status'],fill='white')
    draw.text((10,742),'48 arcsec per panel; Gaia position at center. North up, east left.',fill='#111111')
    sheet.save(ROOT/'optical_unwise.jpg',quality=92)
    (ROOT/'viewer_queries.json').write_text(json.dumps([r for r,_ in results],indent=2)+'\n')

if __name__ == '__main__':
    main()
