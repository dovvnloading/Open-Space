"""Array algebra regressions for the manual workflow, using synthetic data only."""
import numpy as np

from array_search import Cell, Contrast, jackknife, ladder, patterns, sky_grid


def main():
    rng = np.random.default_rng(4821)
    n = 14
    xyz = rng.normal(size=(n,3))
    xyz /= np.linalg.norm(xyz,axis=1)[:,None]
    weight = rng.uniform(0.5,2.,n)
    common_modes = rng.normal(size=(n,4))*0.1
    c = Cell(list(map(str,range(n))),xyz,weight,rng.normal(size=n),
             common_modes @ common_modes.T,np.ones(n),rng.uniform(.4,1.6,n),
             np.zeros((n,5)),np.zeros((n,5)),np.zeros((n,5)))
    for family in ("clock","ephemeris","memory"):
        source = sky_grid()[3] if family == "memory" else None
        g, nuisance = patterns(xyz,family,source)
        b = rng.normal(size=g.shape[1])
        z_signal = weight*(g @ b + nuisance @ rng.normal(size=nuisance.shape[1]))
        model = Contrast(c,family,source)
        np.testing.assert_allclose(model.estimate(z_signal,c.rss).coefficients,b,
                                   rtol=1e-9,atol=1e-10)
        data = np.stack([c.z,c.z+z_signal])
        variances = np.stack([c.rss,c.rss+0.1])
        _, q = ladder(model,data,variances)
        jack = jackknife(c,family,source,data,variances)
        for j,(z,v) in enumerate(zip(data,variances)):
            np.testing.assert_allclose(q[j],model.estimate(z,v).statistic)
            for i in range(n):
                keep = np.arange(n) != i
                sub = Contrast(c.subset(keep),family,source)
                direct = sub.estimate(z[keep],v[keep]).statistic
                np.testing.assert_allclose(jack[i,j],direct,rtol=1e-8,atol=1e-10)


if __name__ == "__main__":
    main()
