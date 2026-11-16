"""Retrieve official DR1 tables, release metadata and image headers; optionally retrieve full images.

python acquire_lodess.py PATH [--images P128+36 P238+32]
Requires astropy and lxml. All URLs are read from the release page.
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

from astropy.io import fits
from lxml import html

RELEASE = 'https://lofar-surveys.org/lodess.html'
PAPER = 'https://arxiv.org/html/2610.04698v1'


def retrieve(url, destination, header_only=False):
    headers = {'Range': 'bytes=0-23039'} if header_only else {}
    for attempt in range(3):
        path = None
        try:
            request = Request(url, headers=headers)
            with urlopen(request, timeout=160) as response:
                status = response.status
                expected = response.headers.get('Content-Length')
                if header_only and (status != 206 or not response.headers.get('Content-Range', '').startswith('bytes 0-')):
                    raise IOError(f'Invalid byte-range response for {url}')
                if not header_only and status != 200:
                    raise IOError(f'Expected whole file, got HTTP {status} for {url}')
                with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as stream:
                    path = Path(stream.name)
                    shutil.copyfileobj(response, stream)
            if expected and path.stat().st_size != int(expected):
                raise IOError(f'Incomplete transfer for {url}')
            if header_only:
                with path.open('rb') as stream:
                    header = fits.Header.fromfile(stream)
                    assert header['SIMPLE'] and header['NAXIS'] >= 2
                    length = stream.tell()
                with path.open('rb+') as stream:
                    stream.truncate(length)
            elif destination.suffix == '.fits':
                with fits.open(path) as hdus:
                    hdus.verify('exception')
                    for hdu in hdus:
                        i = hdu.fileinfo()
                        if i['datLoc'] + i['datSpan'] > path.stat().st_size:
                            raise IOError(f'Incomplete FITS data for {url}')
            os.replace(path, destination)
            return
        except Exception:
            if path:
                path.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def capture(root, images):
    root.mkdir(parents=True, exist_ok=True)
    for folder in ['headers', 'images']:
        (root/folder).mkdir(exist_ok=True)
    for url, filename in [(RELEASE, 'release.html'), (PAPER, 'cordun.html')]:
        retrieve(url, root/filename)
    release = html.parse(str(root/'release.html'))
    urls = [urljoin(RELEASE, u) for u in release.xpath('//a/@href') if u.endswith('.fits')]
    by_field = {u.split('/')[-2]: u for u in urls if '-MFS-I-image-pb.fits' in u}
    qa = html.parse(str(root/'cordun.html'))
    table = next(t for t in qa.xpath('//table') if t.xpath('.//tr[1]/*[normalize-space(.)="Quality"]'))
    rows = [[' '.join(' '.join(c.itertext()).split()) for c in row.xpath('./td|./th')] for row in table.xpath('.//tr')[1:]]
    if {r[0] for r in rows} != set(by_field) or len(rows) != 40:
        raise ValueError('Field/image listing and expected DR1 QA schema do not match.')
    with (root/'fields.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['field', 'ra_deg', 'dec_deg', 'hours', 'flagged_pct', 'rms_mjy_beam', 'quality'])
        writer.writerows(rows)
    targets = [(u, root/('sources.fits' if u.endswith('srl.fits') else 'gaussians.fits'), False)
               for u in urls if u.endswith(('srl.fits', 'gaul.fits'))]
    targets += [(u, root/'headers'/(field+'.hdr'), True) for field, u in by_field.items()]
    targets += [(by_field[f], root/'images'/(f+'-I.fits'), False) for f in images]
    def get(item):
        url, path, partial = item
        retrieve(url, path, partial)
        return {'file': str(path.relative_to(root)), 'url': url,
                'content': 'primary FITS header' if partial else 'complete FITS',
                'bytes': path.stat().st_size, 'sha256': digest(path)}
    with ThreadPoolExecutor(max_workers=6) as pool:
        records = list(pool.map(get, targets))
    with (root/'acquired.json').open('w') as stream:
        json.dump(records, stream, indent=2)
        stream.write('\n')
    (root/'image_urls.txt').write_text('\n'.join(urls)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--images', nargs='*', default=[])
    args = parser.parse_args()
    capture(args.root, args.images)
