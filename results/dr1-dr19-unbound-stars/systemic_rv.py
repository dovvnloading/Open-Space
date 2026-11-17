"""Measure the systemic-RV change required to reverse the strongest escape classification."""

import numpy as np
from infer import draws, escape, load_inputs
from make_queue import save, number
from screen_desi import cartesian

if __name__ == "__main__":
    stars, gs = load_inputs()
    s = next(s for s in stars if s["source_id"] == "1717410430709663616")
    g = gs[s["source_id"]]
    ans = []
    sky = (number(g["ra"]), number(g["dec"]), 1.0, 0.0, 0.0)
    unit = cartesian(*sky, 1.0)[1] - cartesian(*sky, 0.0)[1]
    solar = cartesian(*sky, 0.0)[1]
    solar_proj = float(unit @ solar)
    for model in ("spec_plx", "spdist_only", "bj_photogeo"):
        _, x, v = draws(g, s, model)
        vp = unit @ v
        vt2 = np.sum(v * v, axis=0) - vp**2
        for label, esc in zip(("McMillan17", "Gala_v2"), escape(x)):
            ok = esc**2 > vt2
            lim = np.sqrt(esc[ok] ** 2 - vt2[ok]) - solar_proj
            out = dict(
                source_id=s["source_id"],
                model=model,
                potential=label,
                solar_projection=solar_proj,
                measured_gsr_los=number(s["rv"]) + solar_proj,
                frac_reachable=float(np.mean(ok)),
            )
            for q, val in zip((16, 50, 84), np.percentile(lim, [16, 50, 84])):
                out[f"max_bound_heliocentric_rv_{q}"] = float(val)
            ans.append(out)
    save("systemic_rv.csv", ans)
