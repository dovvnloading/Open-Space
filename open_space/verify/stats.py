"""Independent estimators, written against the textbook definitions.

These exist so that a claim can be recomputed without the code that produced
it. They use only the standard library, take plain sequences of floats, and
favour the obvious formula over speed.
"""

from __future__ import annotations

import datetime
import math

# ---------------------------------------------------------------- descriptive


def finite(values) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(value)]


def mean(values) -> float:
    data = finite(values)
    if not data:
        raise ValueError("mean of no values")
    return math.fsum(data) / len(data)


def median(values) -> float:
    data = sorted(finite(values))
    if not data:
        raise ValueError("median of no values")
    middle = len(data) // 2
    return data[middle] if len(data) % 2 else (data[middle - 1] + data[middle]) / 2


def variance(values, ddof: int = 1) -> float:
    data = finite(values)
    if len(data) <= ddof:
        raise ValueError("variance needs more values than degrees of freedom removed")
    centre = math.fsum(data) / len(data)
    return math.fsum((value - centre) ** 2 for value in data) / (len(data) - ddof)


def std(values, ddof: int = 1) -> float:
    return math.sqrt(variance(values, ddof))


def rms(values) -> float:
    data = finite(values)
    if not data:
        raise ValueError("rms of no values")
    return math.sqrt(math.fsum(value * value for value in data) / len(data))


def mad(values) -> float:
    """Median absolute deviation from the median (not scaled)."""
    data = finite(values)
    centre = median(data)
    return median([abs(value - centre) for value in data])


def quantile(values, q: float) -> float:
    """Linear interpolation between order statistics, the common 'type 7' definition."""
    data = sorted(finite(values))
    if not data:
        raise ValueError("quantile of no values")
    if not 0 <= q <= 1:
        raise ValueError("quantile must lie in [0, 1]")
    position = q * (len(data) - 1)
    low = math.floor(position)
    high = min(low + 1, len(data) - 1)
    return data[low] + (data[high] - data[low]) * (position - low)


# ---------------------------------------------------------------- association


def _paired(x, y) -> tuple[list[float], list[float]]:
    if len(x) != len(y):
        raise ValueError(f"paired samples differ in length: {len(x)} and {len(y)}")
    pairs = [(float(a), float(b)) for a, b in zip(x, y)
             if a is not None and b is not None and math.isfinite(a) and math.isfinite(b)]
    return [a for a, _ in pairs], [b for _, b in pairs]


def pearson(x, y) -> float:
    xs, ys = _paired(x, y)
    if len(xs) < 2:
        raise ValueError("correlation needs at least two pairs")
    mx, my = math.fsum(xs) / len(xs), math.fsum(ys) / len(ys)
    sxy = math.fsum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = math.fsum((a - mx) ** 2 for a in xs)
    syy = math.fsum((b - my) ** 2 for b in ys)
    if sxx == 0 or syy == 0:
        raise ValueError("correlation is undefined for a constant sample")
    return sxy / math.sqrt(sxx * syy)


def ranks(values) -> list[float]:
    """Average ranks, ties sharing the mean of the positions they occupy."""
    order = sorted(range(len(values)), key=lambda index: values[index])
    out = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        shared = (start + stop) / 2 + 1
        for position in range(start, stop + 1):
            out[order[position]] = shared
        start = stop + 1
    return out


def spearman(x, y) -> float:
    xs, ys = _paired(x, y)
    return pearson(ranks(xs), ranks(ys))


def linear_fit(x, y, weights=None) -> tuple[float, float]:
    """Least-squares slope and intercept, optionally weighted (weights multiply squared residuals)."""
    xs, ys = _paired(x, y)
    if weights is None:
        ws = [1.0] * len(xs)
    else:
        if len(weights) != len(x):
            raise ValueError("weights differ in length from the data")
        ws = [float(w) for a, b, w in zip(x, y, weights)
              if a is not None and b is not None and math.isfinite(a) and math.isfinite(b)]
    total = math.fsum(ws)
    if len(xs) < 2 or total <= 0:
        raise ValueError("a line needs at least two weighted points")
    mx = math.fsum(w * a for w, a in zip(ws, xs)) / total
    my = math.fsum(w * b for w, b in zip(ws, ys)) / total
    sxx = math.fsum(w * (a - mx) ** 2 for w, a in zip(ws, xs))
    if sxx == 0:
        raise ValueError("a line is undefined when every x is equal")
    slope = math.fsum(w * (a - mx) * (b - my) for w, a, b in zip(ws, xs, ys)) / sxx
    return slope, my - slope * mx


def weighted_mean(values, weights) -> float:
    xs, ws = _paired(values, weights)
    total = math.fsum(ws)
    if total == 0:
        raise ValueError("weights sum to zero")
    return math.fsum(a * w for a, w in zip(xs, ws)) / total


def chi_square(observed, model, error) -> float:
    if not len(observed) == len(model) == len(error):
        raise ValueError("chi-square inputs differ in length")
    return math.fsum(((o - m) / e) ** 2 for o, m, e in zip(observed, model, error))


# ---------------------------------------------------------------- distributions


def norm_cdf(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def norm_sf(z: float) -> float:
    """One-sided tail probability of a standard normal variable beyond ``z``."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def norm_isf(p: float) -> float:
    """The z whose one-sided tail probability is ``p``. Solved by bisection on ``norm_sf``."""
    if not 0 < p < 1:
        raise ValueError("tail probability must lie strictly between 0 and 1")
    low, high = -40.0, 40.0
    for _ in range(200):
        middle = (low + high) / 2
        if norm_sf(middle) > p:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def _lower_gamma_series(a: float, x: float) -> float:
    term = total = 1.0 / a
    for n in range(1, 10_000):
        term *= x / (a + n)
        total += term
        if abs(term) < abs(total) * 1e-16:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _upper_gamma_fraction(a: float, x: float) -> float:
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for n in range(1, 10_000):
        an = -n * (n - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def gamma_q(a: float, x: float) -> float:
    """Regularised upper incomplete gamma function Q(a, x)."""
    if a <= 0 or x < 0:
        raise ValueError("gamma_q needs a > 0 and x >= 0")
    if x == 0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _lower_gamma_series(a, x)
    return _upper_gamma_fraction(a, x)


def chi2_sf(x: float, k: float) -> float:
    """Probability that a chi-square variable with ``k`` degrees of freedom exceeds ``x``."""
    return gamma_q(k / 2.0, x / 2.0)


def binom_pmf(k: int, n: int, p: float) -> float:
    if not 0 <= k <= n:
        return 0.0
    if p in (0.0, 1.0):
        return 1.0 if k == (n if p == 1.0 else 0) else 0.0
    log_term = (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                + k * math.log(p) + (n - k) * math.log1p(-p))
    return math.exp(log_term)


def binom_sf(k: int, n: int, p: float) -> float:
    """P[X >= k] for X binomial(n, p)."""
    return min(1.0, math.fsum(binom_pmf(i, n, p) for i in range(max(k, 0), n + 1)))


def binom_cdf(k: int, n: int, p: float) -> float:
    """P[X <= k] for X binomial(n, p)."""
    return min(1.0, math.fsum(binom_pmf(i, n, p) for i in range(0, min(k, n) + 1)))


def poisson_pmf(k: int, mu: float) -> float:
    if k < 0:
        return 0.0
    if mu == 0:
        return 1.0 if k == 0 else 0.0
    return math.exp(k * math.log(mu) - mu - math.lgamma(k + 1))


def poisson_sf(k: int, mu: float) -> float:
    """P[X >= k] for X Poisson(mu)."""
    if k <= 0:
        return 1.0
    return 1.0 - gamma_q(float(k), mu)      # P[X >= k] is the regularised lower incomplete gamma P(k, mu)


def poisson_cdf(k: int, mu: float) -> float:
    """P[X <= k] for X Poisson(mu)."""
    return gamma_q(float(k + 1), mu) if k >= 0 else 0.0


# ---------------------------------------------------------------- periodicity


def lomb_scargle_power(times, values, period: float) -> float:
    """Normalised Lomb-Scargle power at one period (Scargle 1982), mean removed."""
    ts, ys = _paired(times, values)
    if len(ts) < 3:
        raise ValueError("a periodogram needs at least three points")
    centre = math.fsum(ys) / len(ys)
    spread = math.fsum((y - centre) ** 2 for y in ys) / (len(ys) - 1)
    if spread == 0:
        raise ValueError("a periodogram is undefined for a constant series")
    omega = 2.0 * math.pi / period
    tau = math.atan2(math.fsum(math.sin(2 * omega * t) for t in ts),
                     math.fsum(math.cos(2 * omega * t) for t in ts)) / (2 * omega)
    cs = sn = cc = ss = 0.0
    for t, y in zip(ts, ys):
        cosine, sine = math.cos(omega * (t - tau)), math.sin(omega * (t - tau))
        cs += (y - centre) * cosine
        sn += (y - centre) * sine
        cc += cosine * cosine
        ss += sine * sine
    return (cs * cs / cc + sn * sn / ss) / (2.0 * spread)


def lomb_scargle_peak(times, values, period_min: float, period_max: float, samples: int = 2000) -> float:
    """Period of the highest Lomb-Scargle power between two periods.

    The interval is scanned uniformly in frequency, then the best sample is refined by
    golden-section search between its neighbours.
    """
    if not 0 < period_min < period_max:
        raise ValueError("need 0 < period_min < period_max")
    samples = max(int(samples), 8)
    f_low, f_high = 1.0 / period_max, 1.0 / period_min
    step = (f_high - f_low) / (samples - 1)
    best_index, best_power = 0, -1.0
    for index in range(samples):
        power = lomb_scargle_power(times, values, 1.0 / (f_low + index * step))
        if power > best_power:
            best_index, best_power = index, power
    left = f_low + max(best_index - 1, 0) * step
    right = f_low + min(best_index + 1, samples - 1) * step
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = left, right
    c, d = b - ratio * (b - a), a + ratio * (b - a)
    pc = lomb_scargle_power(times, values, 1.0 / c)
    pd = lomb_scargle_power(times, values, 1.0 / d)
    for _ in range(60):
        if pc > pd:
            b, d, pd = d, c, pc
            c = b - ratio * (b - a)
            pc = lomb_scargle_power(times, values, 1.0 / c)
        else:
            a, c, pc = c, d, pd
            d = a + ratio * (b - a)
            pd = lomb_scargle_power(times, values, 1.0 / d)
        if abs(b - a) <= 1e-13 * max(abs(a), abs(b)):
            break
    return 2.0 / (a + b)


# ---------------------------------------------------------------- sky and time


def angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Great-circle separation in degrees (Vincenty form, stable at every distance)."""
    r1, d1, r2, d2 = (math.radians(value) for value in (ra1, dec1, ra2, dec2))
    delta = r2 - r1
    numerator = math.hypot(math.cos(d2) * math.sin(delta),
                           math.cos(d1) * math.sin(d2) - math.sin(d1) * math.cos(d2) * math.cos(delta))
    denominator = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(delta)
    return math.degrees(math.atan2(numerator, denominator))


def unit_vector(ra: float, dec: float) -> tuple[float, float, float]:
    r, d = math.radians(ra), math.radians(dec)
    return math.cos(d) * math.cos(r), math.cos(d) * math.sin(r), math.sin(d)


# Dates on which TAI-UTC increased by one second, and the offset from that date on (IERS Bulletin C).
LEAP_SECONDS = (
    ("1972-01-01", 10), ("1972-07-01", 11), ("1973-01-01", 12), ("1974-01-01", 13), ("1975-01-01", 14),
    ("1976-01-01", 15), ("1977-01-01", 16), ("1978-01-01", 17), ("1979-01-01", 18), ("1980-01-01", 19),
    ("1981-07-01", 20), ("1982-07-01", 21), ("1983-07-01", 22), ("1985-07-01", 23), ("1988-01-01", 24),
    ("1990-01-01", 25), ("1991-01-01", 26), ("1992-07-01", 27), ("1993-07-01", 28), ("1994-07-01", 29),
    ("1996-01-01", 30), ("1997-07-01", 31), ("1999-01-01", 32), ("2006-01-01", 33), ("2009-01-01", 34),
    ("2012-07-01", 35), ("2015-07-01", 36), ("2017-01-01", 37),
)
TT_MINUS_TAI = 32.184


def _date(iso: str) -> datetime.datetime:
    text = iso.strip().replace("Z", "").replace("T", " ")
    for pattern in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.datetime.strptime(text, pattern)
        except ValueError:
            continue
    raise ValueError(f"not an ISO date: {iso!r}")


def tai_minus_utc(iso: str) -> int:
    """Whole seconds between TAI and UTC on a UTC date from 1972 onward."""
    day = _date(iso).strftime("%Y-%m-%d")
    offset = None
    for start, seconds in LEAP_SECONDS:
        if day >= start:
            offset = seconds
    if offset is None:
        raise ValueError("leap seconds are defined from 1972-01-01")
    return offset


def tt_minus_utc(iso: str) -> float:
    """Seconds by which Terrestrial Time runs ahead of UTC on a UTC date."""
    return tai_minus_utc(iso) + TT_MINUS_TAI


def julian_date(iso: str) -> float:
    """Julian date of a calendar instant, in whatever time scale the string is given in."""
    moment = _date(iso)
    days = moment.toordinal() + 1721424.5
    return days + (moment.hour * 3600 + moment.minute * 60 + moment.second + moment.microsecond / 1e6) / 86400.0


def modified_julian_date(iso: str) -> float:
    return julian_date(iso) - 2400000.5
