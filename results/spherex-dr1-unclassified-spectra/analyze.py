"""Rank measured, aperture-extracted DR1 Fnu spectra and measure spectral-level recovery.

Run OPENBLAS_NUM_THREADS=1 python analyze.py. Thresholds are search heuristics.
No detector/time-domain injections or formal significances are implied.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import lstsq
from scipy.ndimage import median_filter
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
WIDTHS = (1, 3)
THRESHOLD = 5.0
# Approximate rest centers and intrinsic FWHM in microns; a local passband
# broadening is added in quadrature. These are independent nuisance profiles.
# Ice features: McClure+ 2023 (Nature Astronomy 7,431), Sturm+ 2023 A&A 679,A138.
FEATURES = [
    ('H2O_ice', 3.05, .30, -1), ('NH3_ice', 2.93, .10, -1),
    ('CH4_ice', 3.32, .04, -1), ('CO2_ice', 4.27, .040, -1),
    ('CO_ice', 4.67, .022, -1), ('CO_gas_band', 4.67, .15, -1),
    ('CO_overtone', 2.34, .095, -1), ('H2O_vapor_1', 1.40, .15, -1),
    ('H2O_vapor_2', 1.90, .22, -1), ('CH3OH', 3.53, .035, -1),
    ('PAH', 3.29, .045, +1), ('H_Pag', 1.094, 0, +1),
    ('H_Pab', 1.282, 0, +1), ('H_Paa', 1.8756, 0, +1),
    ('H_Brg', 2.1661, 0, +1), ('H_Brb', 2.6259, 0, +1),
    ('H_Bra', 4.0523, 0, +1), ('H_Pfb', 4.6538, 0, +1),
]
TRIAL_AMPLITUDES = (0.10, 0.25, 0.50, 1.0)


def fits_bank(names, library, grid, wave, lo, hi, teff, redshift):
    """Top-hat passband sampling; stellar and galaxy shapes in relative Fnu."""
    sample = lo[:, None] + (hi - lo)[:, None] * (np.arange(11) + .5) / 11
    ext = wave**-1.7
    temperatures = np.array(sorted({int(n.split('_')[1]) for n in names if n.startswith('CK04_')}))
    selected_t = temperatures
    if np.isfinite(teff):
        selected_t = temperatures[np.argsort(abs(temperatures - teff))[:2]]
    bank, labels = [], []
    z = redshift if np.isfinite(redshift) else 0
    for j, name in enumerate(names):
        stellar = name.startswith('CK04_')
        if stellar and int(name.split('_')[1]) not in selected_t:
            continue
        rest = sample if stellar else sample / (1.0 + z)
        if np.min(rest) < grid[0] or np.max(rest) > grid[-1]:
            continue
        shape = np.interp(rest, grid, library[:, j]).mean(axis=1)
        for av in ((0, .5, 1.5, 4.) if stellar else (0,)):
            reddened = shape * 10**(-.4 * av * .35 * ext)
            labels.append(f'{name}/Av{av:g}')
            bank.append(reddened)
    # Warm dust continuum for embedded sources; the signed ice profiles below
    # provide the absorption freedom. Not a radiative-transfer YSO library.
    for temperature in (450, 800, 1200, 1800):
        shape = wave**-3 / np.expm1(14387.77 / (wave * temperature))
        bank.append(shape / np.median(shape))
        labels.append(f'warm_dust_{temperature}')
    return np.array(bank), labels


def local_filters(ch, detector, wave, error):
    """Whiten line profiles after projecting out local offset + slope per detector."""
    filters, meta = [], []
    for det in np.unique(detector):
        band = np.flatnonzero(detector == det)
        for i in band:
            for width in WIDTHS:
                # Include contiguous support on both sides; grid gaps and seams
                # do not supply sidebands. Wider filters require wider windows.
                radius = 4 if width == 1 else 6
                ix = band[abs(ch[band] - ch[i]) <= radius]
                if min(ch[i] - ch[ix].min(), ch[ix].max() - ch[i]) < radius:
                    continue
                if len(ix) != 2 * radius + 1:
                    continue
                dw = np.median(np.diff(wave[ix]))
                if dw <= 0:
                    continue
                x = (wave[ix] - wave[i]) / dw
                g = np.exp(-4 * np.log(2) * (x / width)**2) / error[ix]
                nuisance = np.column_stack([np.ones(len(ix)), x]) / error[ix, None]
                g -= nuisance @ lstsq(nuisance, g)[0]
                weight = np.zeros(len(ch))
                weight[ix] = g / np.linalg.norm(g)
                filters.append(weight)
                meta.append((i, width, dw))
    return np.array(filters), meta


class Spectrum:
    def __init__(self, tab, names, library, grid, two_band):
        tab = tab.sort_values('channel')
        if tab.channel.duplicated().any() or not tab.channel.between(1, 102).all():
            raise ValueError('Require one-based unique DR1 channels for each source')
        number = lambda key: pd.to_numeric(tab[key], errors='coerce').to_numpy(float)
        flux, sigma = number('flux'), number('sigma')
        ok = np.isfinite(flux) & np.isfinite(sigma) & (sigma > 0)
        ok &= np.isfinite(number('wave_um')) & ~np.isin(number('channel'), [51, 52])
        ok &= (number('nhits') > 0) & (number('flagged') == 0)
        if two_band:
            ok &= np.isin(number('detector'), [2, 5])
        self.ch = number('channel')[ok].astype(int)
        self.det = (self.ch - 1) // 17 + 1
        self.w = number('wave_um')[ok]
        self.y, self.sigma = flux[ok], sigma[ok]
        self.reason = None
        if (len(self.ch) < (15 if two_band else 35) or
                len(np.unique(self.det)) < (2 if two_band else 3)):
            self.reason = 'coverage'
            return
        if np.count_nonzero(self.y / self.sigma >= 3) < 8:
            self.reason = 'source S/N'
            return
        # Freeze the noise scale for injection paired controls. A 2.5% floor
        # handles small proportional mismatches; it is not a calibration claim.
        smooth = median_filter(self.y, size=5, mode='nearest')
        self.e = np.hypot(self.sigma, .025 * abs(smooth))
        self.filters, self.meta = local_filters(self.ch, self.det, self.w, self.e)
        if not self.meta:
            self.reason = 'no contiguous sidebands'
            return
        self.lo = number('wave_lo_um')[ok]
        self.hi = number('wave_hi_um')[ok]
        if np.any(~np.isfinite(self.lo + self.hi) | (self.hi <= self.lo)):
            raise ValueError('Missing/invalid passband widths on measured channels')
        def metadata(name):
            values = pd.to_numeric(tab.get(name, pd.Series(dtype=float)), errors='coerce').dropna().unique()
            if len(values) > 1:
                raise ValueError(f'Inconsistent {name} for a source')
            return values[0] if len(values) else np.nan
        self.z, self.teff = metadata('redshift'), metadata('teff')
        self.bank, self.labels = fits_bank(names, library, grid, self.w, self.lo, self.hi, self.teff, self.z)
        self.bank /= np.sqrt(np.sum((self.bank / self.e)**2, axis=1))[:, None]
        self.bw = self.bank / self.e
        self.profiles = {}
        for z in {0., self.z if np.isfinite(self.z) else 0.}:
            cols, ids = [], []
            for name, rest, fwhm, sign in FEATURES:
                center = rest * (1+z)
                width = np.hypot(fwhm * (1+z), self.hi - self.lo)
                col = sign * np.exp(-4*np.log(2)*((self.w-center)/width)**2) / self.e
                if np.max(abs(col) * self.e) > .15:
                    cols.append(col / np.linalg.norm(col))
                    ids.append(name)
            matrix = np.column_stack(cols) if cols else np.empty((len(self.w), 0))
            reg = np.vstack([matrix, np.eye(len(cols)) * .15])
            self.profiles[z] = matrix, reg, ids

    def fit(self, y):
        wy = y / self.e
        amps = np.maximum(self.bw @ wy, 0)
        # Robust template/normalization search. Each shape retains one amplitude.
        for _ in range(2):
            res = wy - self.bw * amps[:, None]
            weights = np.minimum(1, 2.5 / np.maximum(abs(res), 1.e-12))
            amps = np.maximum(0, np.sum(weights*self.bw*wy, axis=1) / np.sum(weights*self.bw**2, axis=1))
        dist = abs(wy - self.bw * amps[:, None])
        loss = np.where(dist <= 2.5, .5*dist**2, 2.5*(dist - 1.25))
        j = np.argmin(loss.sum(axis=1))
        continuum = self.bank[j] * amps[j]
        r = (y - continuum) / self.e
        z = self.z if self.labels[j].startswith('SWIRE_') and np.isfinite(self.z) else 0.
        profile, reg, ids = self.profiles[z]
        if ids:
            a, _ = nnls(reg, np.r_[r, np.zeros(len(ids))])
            r = r - profile @ a
        scores = self.filters @ r
        return continuum, r, scores, self.labels[j], z


def rank_peers(spectra, min_overlap):
    """Mean distance to the three nearest other measured shapes, no template labels."""
    by_id = {}
    for key, spec in spectra.items():
        positive = spec.y[spec.y > 0]
        shape = np.full(103, np.nan)
        shape[spec.ch] = np.arcsinh(spec.y / (.1 * np.median(positive)))
        by_id[key] = shape
    out = {}
    for key, shape in by_id.items():
        values = []
        for other, peer in by_id.items():
            if key == other:
                continue
            diff = abs(shape - peer)
            if np.isfinite(diff).sum() >= min_overlap:
                med, high = np.nanpercentile(diff, [50, 90])
                values.append(((med + high)*.5, other))
        near = sorted(values)[:3]
        out[key] = (np.mean([d for d, _ in near]) if len(near) == 3 else np.nan,
                    '|'.join(k for _, k in near))
    return out


def inject(key, spec, baseline, rng, repeats):
    continuum, _, sham, _, _ = baseline
    rows = []
    for det in np.unique(spec.det):
        for width in WIDTHS:
            possible = [k for k, (i, wid, dw) in enumerate(spec.meta)
                        if wid == width and spec.det[i] == det and continuum[i] > 0]
            if not possible:
                continue
            for k in rng.choice(possible, size=min(repeats, len(possible)), replace=False):
                i, _, dw = spec.meta[k]
                g = np.exp(-4*np.log(2) * ((spec.w - spec.w[i]) / (width*dw))**2)
                for sign in (-1, 1):
                    base = int(np.argmax(abs(sham)))
                    for fraction in TRIAL_AMPLITUDES:
                        _, _, changed, _, _ = spec.fit(spec.y + sign*fraction*continuum[i]*g)
                        winner = int(np.argmax(abs(changed)))
                        def localized(scores, winner):
                            idx = spec.meta[winner][0]
                            return (spec.det[idx] == det and abs(spec.ch[idx]-spec.ch[i]) <= max(1, width//2)
                                    and sign*scores[winner] >= THRESHOLD)
                        coincident = localized(sham, base)
                        rows.append((key, det, spec.ch[i], spec.w[i], width, sign, fraction,
                                     continuum[i], spec.sigma[i], spec.ch[spec.meta[winner][0]],
                                     changed[winner], int(coincident), int(localized(changed, winner) and not coincident)))
    return rows


def write_frame(path, rows, columns):
    f = pd.DataFrame(rows, columns=columns)
    f.to_csv(path, index=False, float_format='%.7g')
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=HERE / 'spectra.csv')
    ap.add_argument('--templates', type=Path, default=HERE / 'templates/relative_fnu.csv')
    ap.add_argument('--redshifts', type=Path, help='Optional secure source_id,redshift CSV')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--two-band', action='store_true', help='Secondary NEP D2/D5 only; minimum 15 channels in both')
    ap.add_argument('--seed', type=int, default=663931)
    ap.add_argument('--centers-per-detector', type=int, default=2)
    args = ap.parse_args()
    out = args.out
    out.mkdir(exist_ok=True, parents=True)
    data, templates = pd.read_csv(args.input), pd.read_csv(args.templates)
    if args.redshifts:
        extra = pd.read_csv(args.redshifts)
        if 'rvz_redshift' in extra:
            extra = extra.loc[(extra.n_within_3arcsec == 1) & (extra.rvz_type == 'z')].rename(columns={'rvz_redshift':'redshift'})
        extra = extra[['source_id', 'redshift']]
        data = data.drop(columns='redshift', errors='ignore').merge(extra, how='left', on='source_id', validate='m:1')
    all_spectra = {key: Spectrum(t, list(templates.columns[1:]), templates.iloc[:, 1:].to_numpy(),
                                templates.wave_um.to_numpy(), args.two_band)
                   for key, t in data.groupby('source_id', sort=True)}
    write_frame(out/'eligibility.csv', [(key, len(s.ch), len(np.unique(s.det)), s.reason or 'eligible')
                                         for key,s in all_spectra.items()],
                ['source_id', 'n_channels', 'n_detectors', 'status'])
    spectra = {k: s for k,s in all_spectra.items() if s.reason is None}
    peers = rank_peers(spectra, 15 if args.two_band else 35)
    rng = np.random.default_rng(args.seed)
    stats, channels, filters, trials = [], [], [], []
    for key, s in spectra.items():
        baseline = s.fit(s.y)
        continuum, r, scores, label, z = baseline
        win = np.argmax(abs(scores))
        i, wid, _ = s.meta[win]
        model = s.y - r*s.e
        stats.append((key, len(s.ch), len(s.meta), label, s.teff, z, s.ch[i], s.w[i], wid,
                      scores[win], abs(scores[win]), np.median(abs(r)), *peers[key]))
        channels.extend((key, c, w, y, e, m, v) for c,w,y,e,m,v in zip(s.ch,s.w,s.y,s.e,model,r))
        for score, (idx, width, _) in zip(scores, s.meta):
            filters.append((key, s.det[idx], s.ch[idx], s.w[idx], width, score))
        trials.extend(inject(key, s, baseline, rng, args.centers_per_detector))
    stats = pd.DataFrame(stats, columns=['source_id','n_channels','n_filters','template','gaia_teff','model_z','peak_channel',
                                        'peak_wave_um','peak_fwhm_ch','peak_signed','score','median_abs_residual','peer_score','peers'])
    stats['template_rank'] = stats.score.rank(ascending=False, method='min')
    stats['peer_rank'] = stats.peer_score.rank(ascending=False, method='min')
    stats.sort_values('template_rank').to_csv(out/'ranks.csv', index=False, float_format='%.7g')
    write_frame(out/'residuals.csv', channels, ['source_id','channel','wave_um','flux','effective_scale','model','residual_over_scale'])
    filt = write_frame(out/'filter_scores.csv', filters, ['source_id','detector','channel','wave_um','fwhm_ch','score'])
    # Per-center, per-width signed recurrence in the observed frame. Count all
    # eligible other source spectra with testable sidebands, not just winners.
    if not filt.empty:
        f = filt.copy()
        f['ge3'], f['le_neg3'] = f.score >= 3, f.score <= -3
        counts = f.groupby(['channel','wave_um','fwhm_ch']).agg(n_tested=('score','size'),
                        n_emission_ge3=('ge3','sum'), n_absorption_le_neg3=('le_neg3','sum'),
                        median_score=('score','median')).reset_index()
        counts.to_csv(out/'recurrence.csv', index=False, float_format='%.7g')
    t = write_frame(out/'injections.csv', trials, ['source_id','detector','channel','wave_um','fwhm_ch','sign',
                          'amplitude_frac','continuum','sigma','winner_channel','winner_score','sham_coincidence','recovered'])
    if not t.empty:
        for name, keys in [('by_shape',['fwhm_ch','sign','amplitude_frac']),
                           ('by_detector',['detector','fwhm_ch','sign','amplitude_frac'])]:
            grouped = t.groupby(keys).agg(n=('recovered','size'), recovered=('recovered','sum'),
                coincident_shams=('sham_coincidence','sum'), recovery_fraction=('recovered','mean'),
                n_sources=('source_id','nunique')).reset_index()
            grouped.to_csv(out/f'recovery_{name}.csv', index=False, float_format='%.7g')
    input_paths = [args.input, args.templates, Path(__file__)]
    if args.redshifts:
        input_paths.append(args.redshifts)
    receipt = {'input_sha256': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths},
        'two_band':args.two_band,'sources':len(all_spectra),'eligible':len(spectra),
        'n_injections':len(trials),'seed':args.seed,'centers_per_detector':args.centers_per_detector,
        'screening_threshold':THRESHOLD, 'screen_positive':int((stats.score >= THRESHOLD).sum())}
    (out/'analysis.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
