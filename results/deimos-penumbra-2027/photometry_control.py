"""Measure Deimos in the archived Curiosity sol-351 occultation movie."""
import argparse
import csv
import hashlib
import io
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from scipy.ndimage import center_of_mass, find_objects, gaussian_filter, label

BASE = 'https://planetarydata.jpl.nasa.gov/img/data/msl/MSLMST_0004/DATA/RDR/SURFACE/0351/'
# All full-size M products from sequence 001423 on volume MSLMST_0004.
FRAMES = [1172, *range(1191,1199), *range(1205,1213), 1231, 1251]
HEADER = 'product sha256 center_y_px center_x_px phobos_box_distance_px sky_level_RGB_mean_DN flux_r6 flux_r8 flux_r10'.split()


def measure(frame):
    product = f'0351MR001423000030{frame}M01_DRCL'
    with urlopen(BASE+product+'.IMG', timeout=60) as response:
        payload = response.read()
    image = np.frombuffer(payload,dtype=np.uint8).reshape(3,512,512).astype(float).mean(axis=0)
    # Linear track from the well-separated endpoint frames, refined within 6 camera pixels.
    guess_y = 243.5 - 21*(frame-1172)/79
    y,x = np.indices(image.shape)
    radius = np.hypot(y-guess_y,x-222.0)
    background = np.median(image[(radius>11)&(radius<18)])
    cy,cx = center_of_mass(np.maximum(0,image-background)*(radius<6))
    radius = np.hypot(y-cy,x-cx)
    sky = np.median(image[(radius>12)&(radius<18)])
    fluxes = [(image[radius<a]-sky).sum() for a in (6,8,10)]

    # Nearest distance to the box of the tallest bright component (Phobos).
    # Mask the artificial frame edges in the calibrated 8-bit images.
    trimmed = image.copy()
    trimmed[:10] = trimmed[-10:] = 0
    trimmed[:,:10] = trimmed[:,-10:] = 0
    features,_ = label(gaussian_filter(trimmed,1)>50)
    ys,xs = max(find_objects(features), key=lambda obj: obj[0].stop-obj[0].start)
    distance = np.hypot(max(0,ys.start-cy,cy-ys.stop), max(0,xs.start-cx,cx-xs.stop))
    return [product,hashlib.sha256(payload).hexdigest(),cy,cx,distance,sky,*fluxes]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path,help='CSV destination')
    args = p.parse_args()
    with ThreadPoolExecutor(4) as pool:
        rows=list(pool.map(measure,FRAMES))
    text=io.StringIO()
    writer=csv.writer(text,lineterminator='\n')
    writer.writerow(HEADER)
    for row in rows:
        writer.writerow([f'{x:.6f}' if isinstance(x,float) else x for x in row])
    reference=Path(__file__).resolve().with_name('mastcam-control.csv')
    if reference.exists():
        old=list(csv.DictReader(reference.read_text().splitlines()))
        if [(r['product'],r['sha256']) for r in old] != [(r[0],r[1]) for r in rows]:
            raise ValueError('An archived image differs from the pinned sequence')
    args.output.write_text(text.getvalue())
    clean=np.array([r[-3:] for r in rows if r[4]>=20])
    print(f'{len(clean)}/{len(rows)} clear frames. Sample RMS/mean for apertures 6, 8, 10 px: '
          + ', '.join(f'{100*f:.2f}%' for f in clean.std(axis=0,ddof=1)/clean.mean(axis=0)))


if __name__=='__main__':
    main()
