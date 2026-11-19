import sys, numpy as np
from astropy.io import fits
src, fbk = sys.argv[1], sys.argv[2]
h = fits.open(src); lc = h[1]; d = lc.data; hd = lc.header
print('energy band', hd.get('CHANMIN'), hd.get('CHANMAX'), 'TSTART', hd.get('TSTART'), 'ONTIME', h[0].header.get('ONTIME') )
t = d['TIME']; t0 = t[0]
r = np.nan_to_num(d['RATE']); bk = np.nan_to_num(d['BACKV']); fe = d['FRACEXP']
b = np.floor((t - t0) / 300).astype(int); n = len(np.unique(b))
def binned(x, w=None):
    good = fe > 0.5
    s = np.bincount(b[good], weights=x[good], minlength=b.max()+1); c = np.bincount(b[good], minlength=b.max()+1)
    return np.where(c > 0, s / np.maximum(c, 1), np.nan), c
S, c = binned(r); B, _ = binned(bk)
fb = fits.open(fbk)[1]; ft = fb.data['TIME']; fr = np.nan_to_num(fb.data[fb.columns.names[1]])
fbb = np.floor((ft - t0) / 300).astype(int); ok = (fbb >= 0) & (fbb <= b.max())
F = np.bincount(fbb[ok], weights=fr[ok], minlength=b.max()+1) / np.maximum(np.bincount(fbb[ok], minlength=b.max()+1), 1)
hrs = np.arange(len(S)) * 300 / 3600
m = np.isfinite(S); med = np.nanmedian(S); mad = np.nanmedian(np.abs(S[m] - med)) * 1.4826
print(f'duration {hrs[-1]:.2f} h; net rate median {med:.4f} c/s, robust sigma {mad:.4f}; particle fbk median {np.median(F):.3f}')
print('  hr   net     bkg(scaled)  fbk(particle)')
for i in range(len(S)):
    hi = np.isfinite(S[i]) and S[i] > med + 3.5 * mad
    if hi or i % 4 == 0:
        print(f'{hrs[i]:5.2f} {S[i]:7.4f} {B[i]:8.4f} {F[i]:9.3f}{"  <== high" if hi else ""}')
