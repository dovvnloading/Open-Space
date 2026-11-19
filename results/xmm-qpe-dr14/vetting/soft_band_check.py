import sys, numpy as np
from astropy.io import fits
src, fbk = sys.argv[1], sys.argv[2]
h = fits.open(src); d = h[1].data; t = d['TIME']; t0 = h[1].header['TSTART']
soft = np.nan_to_num(d['RATE1']) + np.nan_to_num(d['RATE2']) + np.nan_to_num(d['RATE3'])
sback = np.nan_to_num(d['BACK1V']) + np.nan_to_num(d['BACK2V']) + np.nan_to_num(d['BACK3V'])
hard = np.nan_to_num(d['RATE4']) + np.nan_to_num(d['RATE5'])
fe = d['FRACEXP']; good = fe > 0.5
b = np.floor((t - t0) / 300).astype(int); nb = b.max() + 1
def binned(x):
    s = np.bincount(b[good], weights=x[good], minlength=nb); c = np.bincount(b[good], minlength=nb)
    return np.where(c > 0, s / np.maximum(c, 1), np.nan), c
S, c = binned(soft); SB, _ = binned(sback); H, _ = binned(hard)
fb = fits.open(fbk)[1]; ft = fb.data['TIME']; fr = np.nan_to_num(fb.data[fb.columns.names[1]])
fbb = np.floor((ft - t0) / 300).astype(int); ok = (fbb >= 0) & (fbb < nb)
F = np.bincount(fbb[ok], weights=fr[ok], minlength=nb) / np.maximum(np.bincount(fbb[ok], minlength=nb), 1)
hrs = (np.arange(nb) + 0.5) * 300 / 3600
quiet = np.isfinite(S) & (F < 2 * np.median(F))
med = np.nanmedian(S[quiet]); sig = np.nanmedian(np.abs(S[quiet] - med)) * 1.4826
print(f'0.2-2 keV net: median {med:.4f} c/s, robust sigma {sig:.4f} (quiet-particle bins only, n={quiet.sum()}); particle median {np.median(F):.2f}')
print('MJD(TT)       hr     soft    softbkg  hard(2-12)  particle  note')
mjd0 = 50814.0 + t0 / 86400
for i in range(nb):
    if not np.isfinite(S[i]): continue
    z = (S[i] - med) / sig
    if z > 3.5:
        note = 'PARTICLE FLARE' if F[i] > 2 * np.median(F) else 'quiet particle'
        print(f'{mjd0 + hrs[i]/24:.5f} {hrs[i]:6.2f} {S[i]:8.4f} {SB[i]:8.4f} {H[i]:9.4f} {F[i]:9.2f}  z={z:4.1f} {note}')
