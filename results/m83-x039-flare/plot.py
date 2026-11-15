"""Render the photon time series in the report using its companion CSV and JSON."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

base = Path(__file__).resolve().parent
photons = np.genfromtxt(base / 'photons.csv', delimiter=',', names=True)
summary = json.loads((base / 'measurements.json').read_text())['xray']
time = photons['seconds_from_chip7_gti_start'] / 1000
energy = photons['energy_eV'] / 1000
windows = summary['windows']
edges = [row['start_s'] / 1000 for row in windows] + [windows[-1]['stop_s'] / 1000]
width = np.diff(edges)
with plt.rc_context({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False}):
    fig, (rates, energies) = plt.subplots(2, 1, figsize=(8.0, 5.4), sharex=True,
        gridspec_kw={'height_ratios': [0.9, 1.35]}, layout='constrained')
    fig.suptitle('M83 X039  ·  Chandra observation 13202', x=0.095, ha='left', fontsize=15, weight='semibold')
    for key, color, label in [('source', '#135575', 'Source aperture'),
                               ('expected_background', '#ac532a', 'Local background, scaled to aperture')]:
        values = np.array([row[key] for row in windows]) / width
        rates.stairs(values, edges, color=color, label=label, linewidth=2)
    rates.set_ylabel('Counts per ks')
    rates.set_ylim(0, 2.6)
    rates.legend(frameon=False, loc='upper left')
    energies.scatter(time, energy, s=24, alpha=.8, color='#135575', linewidth=0, label='Individual photon')
    for i, row in enumerate(windows):
        energies.hlines(row['source_median_eV']/1000, edges[i], edges[i+1], color='#ac532a',
                        linewidth=2, label='Window median' if i == 0 else None)
    energies.set_yscale('log')
    energies.set_yticks([.5, 1, 2, 4, 8], labels=['0.5','1','2','4','8'])
    energies.set_ylim(.45, 8.3)
    energies.set_xlim(0, edges[-1])
    energies.set_ylabel('Photon energy (keV)')
    energies.set_xlabel('Time from start of the good-time interval (ks)')
    energies.legend(frameon=False, loc='upper left')
    fig.savefig(base / 'x039-flare.png', dpi=165)
