"""Render two observed Stokes-I dynamic spectra with the component integrations."""
from pathlib import Path
import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from waterfalls import CONFIG, get_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('.'))
    parser.add_argument('--out', type=Path, default=Path('pulse-spectra.png'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.3), sharey=True, constrained_layout=True)
    for ax, (event, (digest, peaks, _split)) in zip(axes, CONFIG.items()):
        with h5py.File(get_data(args.cache, event, digest)) as h:
            use = h['good_freq'][:] & h['calibration/good_freq'][:] & h['flag'][:].all(axis=1)
            data = h['data'][:][use].astype(float)
            freq = h['index_map/freqs'][:][use]
            gates = h.attrs['pulse_emission_region']
            dt = h.attrs['res_time'] * 1000
        off = np.ones(data.shape[1], dtype=bool)
        for start, end in gates:
            off[max(0,start-10):end+11] = False
        data -= data[:, off].mean(axis=1, keepdims=True)
        data /= data[:, off].std(axis=1, keepdims=True)
        bins = []
        for edge in np.arange(400, 800, 10):
            channels = (freq >= edge) & (freq < edge+10)
            if channels.sum() < 40:
                bins.append(np.full(data.shape[1], np.nan))
                continue
            band = data[channels].mean(axis=0)
            bins.append(band / band[off].std())
        pixels = np.stack(bins)
        # Light (three-sample) smoothing only along the time axis for display.
        pixels = np.stack([np.convolve(row, np.ones(3)/3, mode='same') for row in pixels])
        x = (np.arange(data.shape[1]+1) - gates[0].mean() - .5)*dt
        cmap = plt.get_cmap('magma').copy()
        cmap.set_bad('#a7a9ae')
        im = ax.pcolormesh(x, np.arange(400, 810, 10), pixels, cmap=cmap,
                           norm=Normalize(-1.2, 2.8), shading='flat', rasterized=True)
        for (start, end), peak in zip(gates, peaks):
            midpoint = (x[start]+x[end+1])/2
            ax.plot(midpoint, peak, 'o', ms=6, mfc='none', mec='white', mew=1)
        ax.set_xlim(x[gates[0,0]]-10, x[gates[-1,1]]+10)
        ax.set_xlabel('Time from first component (ms)')
        ax.set_title(event.replace('FRB','FRB '), loc='left', fontweight='bold')
    axes[0].set_ylabel('Observed frequency (MHz)')
    fig.colorbar(im, ax=axes, shrink=.85, label='Normalized intensity')
    fig.savefig(args.out, dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    main()
