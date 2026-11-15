"""Check consecutive component centers in Table 2 of Sheikh et al. (2024)."""
import argparse
import csv
from hashlib import sha256
from html import unescape
from pathlib import Path
import re
from urllib.request import urlopen

SOURCE = 'https://arxiv.org/html/2312.07756v1'
DIGEST = '82b582f87ba7e518d75b3402ab2b8ca106fedca933c8a725374770bfba15e916'


def parse_table(raw):
    for table in re.findall(r'<table\b[^>]*>(.*?)</table>', raw, re.S):
        if 'Center Frequency' in table and 'Bandwidth' in table:
            break
    else:
        raise ValueError('Missing published component table.')
    events = {}
    for row in re.findall(r'<tr\b[^>]*>(.*?)</tr>', table, re.S):
        cells = [unescape(re.sub(r'<[^>]*>', '', c)).strip() for c in
                 re.findall(r'<td\b[^>]*>(.*?)</td>', row, re.S)]
        if len(cells) != 9 or not re.match(r'^\d', cells[7]):
            continue
        if cells[0]:
            index = int(cells[0])
            events[index] = (float(cells[1]), [])
        parts = cells[7].split()
        events[index][1].append((float(parts[0]), float(parts[-1])))
    if len(events) != 35 or sum(len(values) for _, values in events.values()) != 101:
        raise ValueError('Unexpected row count.')
    return events


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paper', type=Path, default=Path('ata-paper.html'))
    parser.add_argument('--out', type=Path, default=Path('ata-triplets.csv'))
    args = parser.parse_args()
    if not args.paper.exists():
        args.paper.write_bytes(urlopen(SOURCE, timeout=60).read())
    raw = args.paper.read_bytes()
    if sha256(raw).hexdigest() != DIGEST:
        raise ValueError('Input differs from the arXiv v1 HTML used here.')
    with args.out.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['ATA_index','MJD','triplet_start','center1_MHz','center2_MHz',
                         'center3_MHz','error1_MHz','error2_MHz','error3_MHz','excursion_MHz'])
        for index, (mjd, peaks) in parse_table(raw.decode()).items():
            for start, (a, b, c) in enumerate(zip(peaks, peaks[1:], peaks[2:]), 1):
                writer.writerow([index, mjd, start, a[0], b[0], c[0], a[1], b[1], c[1],
                                 round(b[0] - max(a[0], c[0]), 2)])


if __name__ == '__main__':
    main()
