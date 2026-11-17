"""Check signed morphology for the three triggers/series needing closer review."""
import argparse
import csv
from pathlib import Path

import numpy as np

from scan import Binning, detect, read_curve, write_csv
from vet import average, read_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    results = []
    cases = {'257556370': (4208.90,4209.75), '394558249': (4223.4,4230.6),
             '412049646': (4220.35,4222.25)}
    for row in read_rows(args.manifest):
        tic = row['tic']
        if tic not in cases:
            continue
        c = read_curve(args.manifest.parent/row['path'], row)
        bins = Binning(c)
        y = bins.bin_flux(bins.flux)
        det = average(bins,c,c.diagnostics['DET_FLUX'])
        small = average(bins,c,c.diagnostics['DET_FLUX_SML'])
        x = average(bins,c,c.diagnostics['SAP_X'])
        sky = average(bins,c,c.diagnostics['SAP_BKG'])
        lo, hi = cases[tic]
        for i in np.flatnonzero((bins.centers > lo) & (bins.centers < hi) & bins.valid):
            results.append(dict(tic=tic,btjd=bins.centers[i],sap=y[i],det=det[i],
                                small=small[i],sap_x=x[i],sap_bkg=sky[i]))
        if tic == '257556370':
            # Keep timestamps, shoulders and slowly varying SAP response. Remove
            # only the localized positive DET pulse above the sector median + 1%.
            use = c.good & np.isfinite(c.time) & np.isfinite(c.flux)
            corrected = c.diagnostics['DET_FLUX'][use]
            cap = np.nanmedian(corrected) + .01
            remove = ((bins.time > 4209.24) & (bins.time < 4209.47) &
                      np.isfinite(corrected) & (corrected > cap))
            trial = bins.flux.copy()
            trial[remove] *= cap / corrected[remove]
            after = detect(bins.bin_flux(trial),bins.centers)
            print(f'{tic}: {int(remove.sum())} positive cadences capped at DET={cap:.6f}; '
                  f'original {len(detect(y,bins.centers))} triggers, {len(after)} after')
    write_csv(args.out,results,list(results[0]))


if __name__ == '__main__':
    main()
