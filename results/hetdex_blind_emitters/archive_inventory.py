"""Read release FITS headers using HTTP ranges; record checksums for the bytes inspected."""
import hashlib
import json
from math import prod
from pathlib import Path
import re
from urllib.request import Request, urlopen

BASE = 'https://web.corral.tacc.utexas.edu/hetdex/HETDEX/'
RELEASES = [
    ('sc1', '3.2', 'catalogs/hetdex_source_catalog_1/', 7448504),
    ('sc2', '1.5', 'pdr/pdr1/hetdex_source_catalog_2/', 19581262),
]


def headers(url):
    offset = 0
    info = []
    while True:
        request = Request(url, headers={'Range': f'bytes={offset}-{offset + 65535}', 'User-Agent': 'Open-Space HETDEX research'})
        with urlopen(request, timeout=50) as response:
            if response.status != 206:
                raise ValueError(f'range returned status {response.status}')
            byte_range = response.headers['Content-Range']
            first, end, total = map(int, re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', byte_range).groups())
            if first != offset:
                raise ValueError(f'incorrect range {byte_range}')
            block = response.read()
        cards = []
        for i in range(0, len(block) - 79, 80):
            card = block[i:i + 80].decode('ascii')
            cards.append(card.rstrip())
            if card[:8].strip() == 'END':
                break
        else:
            raise ValueError('FITS header longer than probe')
        values = {c[:8].strip(): c[10:].split('/')[0].strip().strip("'").strip() for c in cards if c[8:10] == '= '}
        naxis = int(values['NAXIS'])
        dims = [int(values[f'NAXIS{i + 1}']) for i in range(naxis)]
        data_size = (prod(dims) if dims else 0) * abs(int(values['BITPIX'])) // 8 + int(values.get('PCOUNT', 0))
        size = ((len(cards) * 80 + 2879) // 2880) * 2880
        info.append(dict(offset=offset, data_offset=offset + size, cards=cards,
                         header_sha256=hashlib.sha256(block[:size]).hexdigest()))
        offset += size + (data_size + 2879) // 2880 * 2880
        if offset == total:
            return {'url': url, 'bytes': total, 'hdus': info}
        if offset > total or len(info) > 30:
            raise ValueError('unexpected FITS layout')


def main():
    result = {}
    for tag, version, path, record in RELEASES:
        filename = f'hetdex_{tag}_spec_v{version}.fits'
        attempts = []
        for url in [BASE + path + filename, f'https://zenodo.org/records/{record}/files/{filename}?download=1']:
            try:
                entry = headers(url)
            except Exception as exc:
                attempts.append({'url': url, 'error': str(exc)})
            else:
                result[tag] = dict(attempts=attempts, **entry)
                break
        else:
            result[tag] = {'attempts': attempts}
    out = Path(__file__).parent / 'evidence'
    out.mkdir(exist_ok=True)
    (out / 'archive_headers.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
