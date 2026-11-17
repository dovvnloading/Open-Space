"""Leave-one-source-out check of the fixed-observed-channel residual distribution."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis_dir', type=Path)
    args = parser.parse_args()
    directory = args.analysis_dir
    df = pd.read_csv(directory/'filter_scores.csv')
    blocks = []
    for _, group in df.groupby(['channel', 'fwhm_ch']):
        if group.source_id.duplicated().any():
            raise ValueError('Expected one score per source/channel/width')
        if len(group) < 11:
            continue
        group = group.copy()
        scores = group.score.to_numpy()
        peers = np.repeat(scores[None, :], len(group), axis=0)
        np.fill_diagonal(peers, np.nan)
        med = np.nanmedian(peers, axis=1)
        mad = np.nanmedian(abs(peers - med[:,None]), axis=1)
        # The unit floor limits amplification in a fortuitously quiet channel.
        group['adjusted'] = (scores - med) / np.maximum(1., 1.4826*mad)
        group['other_n'] = len(group)-1
        group['other_same_sign_ge3'] = np.sum(peers*np.sign(scores[:,None]) >= 3, axis=1)
        group['other_abs_at_least_target'] = np.sum(abs(peers) >= abs(scores[:,None]), axis=1)
        group['other_median'] = med
        group['other_mad'] = mad
        blocks.append(group)
    if not blocks:
        raise ValueError('No channel/width has ten independent other sources')
    comparisons = pd.concat(blocks)
    ranks = pd.read_csv(directory/'ranks.csv')
    original = ranks[['source_id', 'peak_channel', 'peak_fwhm_ch']].rename(
        columns={'peak_channel':'channel', 'peak_fwhm_ch':'fwhm_ch'})
    original = original.merge(comparisons, on=['source_id','channel','fwhm_ch'], how='left', validate='1:1')
    strongest = comparisons.loc[comparisons.adjusted.abs().groupby(comparisons.source_id).idxmax().to_numpy()].copy()
    strongest['adjusted_rank'] = strongest.adjusted.abs().rank(ascending=False, method='min')
    original.merge(strongest, on='source_id', how='outer', suffixes=('_original','_adjusted'),
                   validate='1:1').to_csv(directory/'peer_channel_check.csv', index=False, float_format='%.7g')


if __name__ == '__main__':
    main()
