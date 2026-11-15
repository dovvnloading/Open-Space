"""Render each GMN track's target thumbnails and adjacent-detection controls locally."""
import argparse
import csv
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from PIL import Image, ImageDraw


def get_mosaic(record, cache):
    path = cache / (record['station'] + '.jpg')
    if path.exists():
        data = path.read_bytes()
    else:
        request = Request(record['detected_mosaic_url'], headers={'User-Agent': 'GMN-camera-review/1.0'})
        with urlopen(request, timeout=40) as response:
            data = response.read()
    if sha256(data).hexdigest() != record['mosaic_sha256']:
        raise ValueError(f"Content changed for {record['station']}")
    if not path.exists():
        path.write_bytes(data)
    return Image.open(BytesIO(data)).convert('L')


def crop(mosaic, index):
    x = 320 * (index % 10)
    y = 30 + 190 * (index // 10)
    return mosaic.crop((x, y, x + 320, y + 180))


def draw_track(records, cache):
    # Each camera gets one row: original block, same-field median of four
    # time neighbours, and a difference to highlight transient streaks.
    sheet = Image.new('L', (960, len(records) * 208 + 20), color=16)
    draw = ImageDraw.Draw(sheet)
    for i, name in enumerate(('Candidate block', 'Median of nearby detected blocks', '4x difference + 32 grey levels')):
        draw.text((3 + 320*i, 3), name, fill=255)
    for i, record in enumerate(records):
        mosaic = get_mosaic(record, cache)
        index = 10 * int(record['mosaic_row_zero_based']) + int(record['mosaic_column_zero_based'])
        original = crop(mosaic, index)
        others = [np.asarray(crop(mosaic, index + delta), dtype=float) for delta in (-2, -1, 1, 2)]
        background = np.median(others, axis=0)
        difference = np.clip(4 * (np.asarray(original) - background) + 32, 0, 255)
        y = 20 + 208 * i
        label = f"{record['station']}  {record['ff_label_utc']} UTC   {record['visual_description']}"
        draw.text((3, y + 4), label, fill=255)
        for j, item in enumerate((original, Image.fromarray(background.astype('uint8')),
                                  Image.fromarray(difference.astype('uint8')))):
            sheet.paste(item, (320*j, y + 20))
    return sheet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    with Path(__file__).with_name('camera-tiles.csv').open(newline='') as stream:
        records = list(csv.DictReader(stream))
    for track in ('1', '2', '3', '4'):
        rows = [r for r in records if r['track'] == track]
        draw_track(rows, args.cache).save(args.output / f'track-{track}.png')


if __name__ == '__main__':
    main()
