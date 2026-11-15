import pandas as pd, numpy as np, itertools

kois = pd.read_csv("koi.csv")
ks = kois[(kois.koi_disposition == "CONFIRMED") & (kois.koi_model_snr > 10)].copy()
ps = []
for k, grp in ks.groupby("kepid"):
    for a, b in itertools.combinations(grp.itertuples(), 2):
        snr = min(x.koi_model_snr * np.sqrt(x.koi_period / 1410) for x in (a, b))
        t = a.koi_time0bk + np.arange(-1, int(1610 / a.koi_period) + 2) * a.koi_period
        t = t[(t > 131) & (t < 1591)]
        u = b.koi_time0bk + np.round((t - b.koi_time0bk) / b.koi_period) * b.koi_period
        ok = abs(t - u) < (a.koi_duration + b.koi_duration) / 48
        if any(ok):
            ps.append(
                [
                    snr,
                    k,
                    a.kepler_name,
                    b.kepler_name,
                    len(t[ok]),
                    round(a.koi_depth),
                    round(b.koi_depth),
                    round(a.koi_impact, 2),
                    round(b.koi_impact, 2),
                ]
            )
df = pd.DataFrame(ps, columns="snr k a b n d1 d2 b1 b2".split()).sort_values("snr", ascending=False)
df.to_csv("pairs.csv", index=False)
print(df.head(32).to_string(index=False))
