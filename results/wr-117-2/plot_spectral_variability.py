"""Plot disjoint wavelength matches; vertical bars give formal flux errors."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

ROOT = Path(__file__).resolve().parent


def main():
    with (ROOT/'six-detector-matches.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), sharex='col',
                             gridspec_kw={'height_ratios': [1.5, 1]}, layout='constrained')
    styles = [('first', '#355c8a', 'o', '2025 Aug–Sep'), ('middle', '#aa4a16', 's', '2026 spring'),
              ('late', '#22765e', '^', '2026 Aug–Oct')]
    for col, detector in enumerate([5, 6]):
        stamps = [r for r in rows if int(r['detector']) == detector]
        a, b = axes[:, col]
        for visit, color, marker, label in styles:
            wave = np.array([float(r[visit+'_um']) for r in stamps])
            flux = np.array([float(r[visit+'_jy']) for r in stamps])
            err = np.array([float(r[visit+'_formal_error_jy']) for r in stamps])
            a.errorbar(wave, flux*1000, err*1000, fmt=marker, markersize=5,
                       color=color, ecolor=color, elinewidth=.9, capsize=2, label=label)
        first = np.array([float(r['first_jy']) for r in stamps])
        a.set_title(f'Detector {detector}   ·   {len(stamps)} triplets', loc='left', fontweight='bold')
        for visit, color, marker, _ in styles[1:]:
            after = np.array([float(r[visit+'_jy']) for r in stamps])
            wave = np.array([float(r[visit+'_um']) for r in stamps])
            error = np.hypot([float(r['first_formal_error_jy']) for r in stamps],
                             [float(r[visit+'_formal_error_jy']) for r in stamps])
            b.errorbar(wave, (after-first)*1000, error*1000, fmt=marker, markersize=5,
                       color=color, ecolor=color, elinewidth=.9, capsize=2)
        b.axhline(0, lw=.8, color='0.4')
        b.set_xlabel('Center wavelength (μm)')
        for axis in (a,b):
            axis.grid(axis='y', alpha=.18)
            axis.spines[['top','right']].set_visible(False)
            axis.yaxis.set_major_locator(MaxNLocator(5))
    axes[0,0].set_ylabel('Aperture flux density (mJy)')
    axes[1,0].set_ylabel('Difference from first visit (mJy)')
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc='outside lower center')
    fig.suptitle('WR 117-2 sightline   |   Matched SPHEREx spectral images',fontsize=15)
    fig.savefig(ROOT/'spectral-variability.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
