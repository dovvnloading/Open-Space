"""Independent, offline equal-area XRT echo checks. Supply the four-row photon .npy."""
import argparse
from pathlib import Path
import numpy as np
from scipy.special import xlogy

BANDS = [.7, .9, 1.1, 1.3]


def increments(data):
    index = np.argsort(data[0])
    visit = np.empty_like(index)
    visit[index] = np.r_[0, np.cumsum(np.diff(data[0, index]) > 500)]
    return visit


def aperture(data, distance, padding=0):
    target = 826.6 * data[0] / distance
    complete = (.7 * target > (50 + padding)**2) & (1.3 * target < (550 - padding)**2)
    return target, complete


def counts(data, distance, subset, center=(0, 0), padding=0):
    target, complete = aperture(data, distance, padding)
    ratio = ((data[1] - center[0])**2 + (data[2] - center[1])**2) / target
    return np.histogram(ratio[complete & subset], BANDS)[0]


def score(count):
    inner,on,outer = count
    total = count.sum()
    if not total or on <= total / 3:
        return 0.
    return np.sqrt(2 * (xlogy(on, 3 * on / total) +
                        xlogy(inner + outer, 1.5 * (inner + outer) / total)))


def late_reference(data, distance, soft, visits):
    target, complete = aperture(data, distance)
    late = soft & (.7 * target > 550**2)
    if not np.any(late):
        return
    squared_radius = np.sort((data[1]**2 + data[2]**2)[late])
    observed = expected = 0.
    print(f'Late fixed-sky control: {late.sum()} photons, {len(np.unique(visits[late]))} snapshots')
    for v in np.unique(visits[soft & complete]):
        use = soft & complete & (visits == v)
        brackets = np.searchsorted(squared_radius, np.asarray(BANDS)[:, None] * target[use])
        template = np.diff(brackets, axis=0).sum(axis=1)
        early = counts(data, distance, use)
        if template.sum() == 0:
            print(v, 'uncovered')
            continue
        pred = template[1] * early.sum() / template.sum()
        print(f'  snapshot {v}: observed {early[1]}, template {pred:.2f}')
        expected += pred
        observed += early[1]
    print(f'Observed {observed:g}, template {expected:.2f} on-annulus photons.')


def validate(path, distance):
    data = np.load(path)
    visits = increments(data)
    soft = (data[3] >= 70) & (data[3] <= 300)
    print('Inner/on/outer, 0.7–3 keV:', counts(data, distance, soft))
    hard = (data[3] >= 301) & (data[3] <= 700)
    print('Inner/on/outer, 3.01–7 keV:', counts(data, distance, hard))
    for label, subset in [('first three', visits < 3), ('later', visits >= 3),
                          ('even', visits % 2 == 0), ('odd', visits % 2 == 1)]:
        distances = np.exp(np.arange(np.log(20), np.log(20000), .015))
        statistics = [score(counts(data, d, soft & subset)) for d in distances]
        best = distances[np.argmax(statistics)]
        print(f'{label}: best {best:.1f} pc, fixed-D counts {counts(data, distance, soft & subset)}')
    octant = np.floor((np.arctan2(data[2],data[1]) + np.pi) * 4 / np.pi).astype(int) % 8
    octants = [counts(data,distance,soft & (octant == i)) for i in range(8)]
    print('Octant on minus mean of sidebands:',[float(b-(a+c)/2) for a,b,c in octants])
    for offset in (20,40,60,90):
        displaced = []
        for angle in np.arange(16) * np.pi / 8:
            center = (offset * np.cos(angle), offset * np.sin(angle))
            displaced.append(score(counts(data, distance, soft, center, offset)))
        centered = score(counts(data, distance, soft, padding=offset))
        print(f'{offset} arcsec: centered {centered:.2f}, max among 16 offset directions {max(displaced):.2f}')
    late_reference(data, distance, soft, visits)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('photons',type=Path)
    parser.add_argument('distance_pc',type=float)
    args = parser.parse_args()
    validate(args.photons, args.distance_pc)
