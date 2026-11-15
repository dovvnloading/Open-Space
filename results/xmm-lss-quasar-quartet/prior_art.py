"""Crossmatch the published pair and multiple-AGN tables against the field center."""
import csv
import hashlib
import io
import math
import urllib.request

RA, DEC = 36.31513, -4.78245
JING = 'https://content.cld.iop.org/journals/0004-637X/1000/2/311/revision1/apjae517at1_mrt.txt'
BIGMAC = 'https://raw.githubusercontent.com/thatastroguy/thebigmac/main/DR1_tabledownloads/BigMAC_{}table_DR1.csv'


def distance(ra, dec):
    half_chord_sq = (math.sin(math.radians((DEC-dec)/2))**2
                     + math.cos(math.radians(DEC))*math.cos(math.radians(dec))
                     * math.sin(math.radians((RA-ra)/2))**2)
    return math.degrees(2*math.asin(math.sqrt(half_chord_sq)))


def sexagesimal(field):
    a, m, s = map(float, field.strip().removesuffix('s').split(':'))
    return (-1 if field.startswith('-') else 1) * (abs(a) + m/60 + s/3600)


if __name__ == '__main__':
    with urllib.request.urlopen(JING, timeout=60) as response:
        raw = response.read()
    assert hashlib.sha256(raw).hexdigest() == '29862aea353860b4208dc1c0ba119e906f4ee41cc2d9617bc761cdc933cf5900'
    rows = [line for line in raw.decode().splitlines() if line.startswith('J')]
    separations = [distance(float(row[start:end]), float(row[start+10:end+10]))
                   for row in rows for start, end in [(67, 76), (301, 310)]]
    print('Jing pairs, nearest component (arcmin):', len(rows), min(separations)*60)
    for kind in ['main', 'triple', 'quad', 'quint']:
        with urllib.request.urlopen(BIGMAC.format(kind), timeout=60) as response:
            rows = list(csv.DictReader(io.StringIO(response.read().decode('utf-8-sig'))))
        separations = []
        for row in rows:
            for i in range(1, 6):
                ra = row.get(f'RA{i}', '')
                dec = row.get(f'Dec{i}', '')
                if ':' in ra and ':' in dec:
                    separations.append(distance(sexagesimal(ra)*15, sexagesimal(dec)))
        print('Big MAC', kind, 'rows, component coordinates, nearest separation (degrees):',
              len(rows), len(separations), min(separations))
