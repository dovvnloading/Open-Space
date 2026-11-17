"""Known Galactic-centre ejecta S5-HVS1. The distance input is independent of an assumed origin."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from orbits import integrate
from infer import summarize
from make_queue import rows, save

if __name__ == "__main__":
    g = rows("positive/gaia/rows_000.csv")[0]
    d = 10**0.936
    # Koposov et al. Table 1 separate from the GC-constrained distance; log10(d/kpc)=0.936 +/- 0.015.
    s = dict(
        source_id=g["source_id"],
        channel="known_S5",
        spec_d=d,
        spec_e=d * np.sqrt(np.expm1((0.015 * np.log(10)) ** 2)),
        rv=1017.0,
        e_rv=2.7,
        zp=0.0,
    )
    result = next(r for r in summarize(s, g) if r["model"] == "spec_only")
    save("positive/recovery.csv", [result])
    save(
        "positive/orbit.csv",
        [integrate(s, g, "spec_only", dt=0.025, n=256, duration=10.0)],
    )
