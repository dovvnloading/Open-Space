"""Transfer each catalog-redshift object's winning feature to the other redshifts."""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('field', type=Path)
    args = parser.parse_args()
    field = args.field
    spectra = pd.read_csv(field/'spectra.csv')
    redshifts = pd.read_csv(field/'catalogs/model_redshifts.csv').set_index('source_id').redshift
    ranks = pd.read_csv(field/'full/ranks.csv').set_index('source_id')
    filters = pd.read_csv(field/'full/filter_scores.csv')
    residuals = pd.read_csv(field/'full/residuals.csv')
    rows = []
    for a, za in redshifts.items():
        if a not in ranks.index:
            continue
        peak = ranks.loc[a]
        for b, zb in redshifts.items():
            if a == b or b not in ranks.index:
                continue
            for frame, expected in [('fixed_observed', peak.peak_wave_um),
                                    ('common_rest', peak.peak_wave_um * (1+zb)/(1+za))]:
                bins = spectra[(spectra.source_id == b) & (spectra.wave_lo_um <= expected) &
                               (spectra.wave_hi_um >= expected)]
                c = bins.loc[(bins.wave_um - expected).abs().idxmin()].channel if len(bins) else None
                result = filters[(filters.source_id == b) & (filters.channel == c) &
                                 (filters.fwhm_ch == peak.peak_fwhm_ch)]
                residual = residuals[(residuals.source_id == b) & (residuals.channel == c)]
                rows.append(dict(from_id=a, from_z=za, from_peak_um=peak.peak_wave_um, from_score=peak.peak_signed,
                                 to_id=b,to_z=zb,prediction=frame,expected_um=expected,
                                 channel=c, n_containing_bins=len(bins), same_width_score=result.score.iloc[0] if len(result) else None,
                                 residual_over_scale=residual.residual_over_scale.iloc[0] if len(residual) else None))
    pd.DataFrame(rows).to_csv(field/'full/redshift_check.csv',index=False,float_format='%.8g')


if __name__ == '__main__':
    main()
