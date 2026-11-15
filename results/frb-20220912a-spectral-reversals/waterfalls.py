"""Measure November 29 pulse colors from the released CHIME Stokes-I arrays."""
import argparse
import csv
from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen

import h5py
import numpy as np

BASE = ('https://cadc-west-01.canfar.net/vault/files/AstroDataCitationDOI/'
        'CISTI.CANFAR/25.0066/data/dynamic_spectra/data/')
CONFIG = {
    'FRB20221129A': ('88d4e3ff4aa4aeed065f025aa59b1884afc9ddc338f68bb73df224ec0420b585',
                     (501.35, 560.95, 503.67), 530.0),
    'FRB20221129B': ('dc2f4b650bb4195b18f3668243450d42c1b5cb7e364b89d2371dc0928d0492b5',
                     (565.91, 706.21, 546.30), 650.0),
}


def get_data(cache, event, digest):
    path = cache / (event + '.h5')
    if not path.exists():
        data = urlopen(BASE + event + '_stokesi_dynamic_spectrum.h5', timeout=120).read()
        if sha256(data).hexdigest() != digest:
            raise ValueError(f'Unexpected download for {event}')
        path.write_bytes(data)
    if sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f'Unexpected input for {event}')
    return path


def measure(path, peaks, split):
    with h5py.File(path) as h:
        # Include only channels usable at every time, also passing both released masks.
        use = h['good_freq'][:] & h['calibration/good_freq'][:] & h['flag'][:].all(axis=1)
        f = h['index_map/freqs'][:][use]
        data = h['data'][:][use].astype(float)
        gates = h.attrs['pulse_emission_region']
        dm = float(h.attrs['dm_incoherent'])
        dt = float(h.attrs['res_time'])
        k = float(h.attrs['dm_const'])
        beam = int(h.attrs['beam_number'])
    off = np.ones(data.shape[1], dtype=bool)
    for start, end in gates:
        off[max(0, start-10):end+11] = False
    data -= data[:, off].mean(axis=1, keepdims=True)
    scale = data[:, off].std(axis=1)
    use = np.isfinite(scale) & (scale > 0)
    f, data = f[use], data[use] / scale[use, None]
    # Integrate all off-pulse gates of each matched duration. Subtract the
    # empirical gate noise from the observed two-band sums to gauge sensitivity.
    bands = np.stack([data[f < split].sum(axis=0), data[f >= split].sum(axis=0)])
    envelopes = []
    for start, end in gates:
        width = end - start + 1
        blank = np.flatnonzero(np.convolve(off.astype(int), np.ones(width, int),
                                          mode='valid') == width)
        noise = np.stack([bands[:, i:i+width].sum(axis=1) for i in blank])
        signal = bands[:, start:end+1].sum(axis=1)
        draws = signal - (noise - noise.mean(axis=0))
        fractions = draws[:, 1] / draws.sum(axis=1)
        envelopes.append([len(blank), *np.round(np.quantile(fractions, [.025, .975]), 5)])
    rows = []
    for trial in dict.fromkeys((round(dm, 6), 219.46, 219.893117, 220.5, 221.5, 222.938288)):
        for number, (peak, (start, end)) in enumerate(zip(peaks, gates), 1):
            # Retain the measured arrival time at the fitted peak as DM varies.
            # Changes in predicted arrival times elsewhere follow cold-plasma dispersion.
            shift = np.rint(k*(trial-dm)*(f**-2 - peak**-2)/dt).astype(int)
            bins = np.arange(start, end+1)[None, :] + shift[:, None]
            if bins.min() < 0 or bins.max() >= data.shape[1]:
                raise ValueError('A gate falls outside the cached data.')
            energy = data[np.arange(len(f))[:, None], bins].sum(axis=1)
            high = f >= split
            # These are fractions of baseline-subtracted, variance-normalized counts.
            # Noise can drive their estimates outside [0,1].
            frac = energy[high].sum() / energy.sum()
            centroid = (f * energy).sum() / energy.sum()
            rows.append([beam, len(f), round(trial, 6), number, start, end, split,
                         round(frac, 5), round(centroid, 3),
                         *(envelopes[number-1] if abs(trial-dm) < 1e-5 else ['', '', ''])])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('.'))
    parser.add_argument('--out', type=Path, default=Path('pulse-colors.csv'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, (digest, peaks, split) in CONFIG.items():
        path = get_data(args.cache, name, digest)
        rows.extend([name, *row] for row in measure(path, peaks, split))
    with args.out.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['event','beam','channels','DM_pc_cm-3','component','first_bin',
                         'last_bin','band_division_MHz','high_band_fraction','centroid_MHz',
                         'off_gates','noise_low','noise_high'])
        writer.writerows(rows)


if __name__ == '__main__':
    main()
