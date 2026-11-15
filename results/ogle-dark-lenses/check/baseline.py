"""Baseline stability of OGLE-2026-BLG-0058 before the event (2017-2024).

Usage: python -I check/baseline.py <phot.dat>
"""
import sys

import numpy as np
from astropy.time import Time
from astropy.timeseries import LombScargle

d = np.loadtxt(sys.argv[1])
t, mag, err, seeing, sky = d.T
years = Time(t, format="jd").decimalyear
print("span %.2f-%.2f, N=%d, median err %.3f" % (years.min(), years.max(), len(t), np.median(err)))
for y in range(int(years.min()), int(years.max()) + 1):
    m = (years >= y) & (years < y + 1)
    if m.sum() == 0:
        continue
    med = np.median(mag[m])
    mad = 1.4826 * np.median(np.abs(mag[m] - med))
    print("%d N=%4d median I=%.3f robust scatter=%.3f median err=%.3f" % (y, m.sum(), med, mad, np.median(err[m])))

base = years < 2024.9
tb, mb, eb = t[base], mag[base], err[base]
coef = np.polyfit(tb - tb.mean(), mb, 1, w=1 / eb)
print("baseline linear trend: %.4f mag per year" % (coef[0] * 365.25))
ls = LombScargle(tb, mb, eb)
freq, power = ls.autopower(minimum_frequency=1 / 1500, maximum_frequency=1 / 1.05, samples_per_peak=10)
top = np.argsort(power)[::-1]
shown = []
for i in top:
    p = 1 / freq[i]
    if all(abs(p - q) / q > 0.05 for q in shown):
        shown.append(p)
        print("peak P=%.2f d power=%.4f FAP=%.2e" % (p, power[i], ls.false_alarm_probability(power[i])))
    if len(shown) == 4:
        break
# outliers relative to baseline
z = (mb - np.median(mb)) / eb
print("baseline points beyond 5 sigma: %d of %d" % ((np.abs(z) > 5).sum(), len(mb)))
