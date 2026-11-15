"""Exercise the illumination, numerical, alternate ephemeris and provenance checks."""
import hashlib
import json
import numpy as np

from replay import INPUT, calculate


def reference_values(cassini=False, alternate=False, limb=.6, ring_scale=1., rays=1024):
    rows = calculate(rays, limb, ring_scale, 5, alternate, cassini)
    lookup = {r[0][5:16]: r[1:] for r in rows}
    return np.array([lookup[t][:2] for t in ["03-01T22:05", "03-02T05:00",
                                             "03-02T16:45", "03-03T04:35", "03-03T10:25"]])


if __name__ == "__main__":
    for meta in INPUT.glob("*_from_*.json"):
        expected = json.loads(meta.read_text())["sha256_file"]
        assert hashlib.sha256(meta.with_suffix(".csv").read_bytes()).hexdigest() == expected
    default = reference_values()
    assert np.all(default[[0, 4], 1] == 1) and default[2, 1] == 0
    for kwargs in [{}, {"cassini": True}, {"alternate": True}, {"limb": 0},
                   {"ring_scale": .5}, {"ring_scale": 2}]:
        f = reference_values(**kwargs)[:, 0]
        assert max(f[0], f[4]) < .15 and min(f[1], f[3]) > .8 and f[2] == 0
    assert np.max(np.abs(reference_values(rays=8192) - default)) < .002
    print("Five input hashes, both measured opacity profiles, PH20, optical-depth "
          "perturbations, a uniform solar disk, and solar-disk quadrature passed.")
