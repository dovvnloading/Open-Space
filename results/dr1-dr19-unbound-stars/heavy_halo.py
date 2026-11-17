"""Published high-mass-halo sensitivity, separated from the two reference Galactic potentials."""

import numpy as np
from galpy.potential.mwpotentials import Irrgang13III
from galpy.potential.interpRZPotential import eval_potential_c
from infer import draws, load_inputs
from make_queue import save

if __name__ == "__main__":
    stars, gs = load_inputs()
    out = []
    pot = Irrgang13III
    for s in stars:
        if s["source_id"] not in ("4444095982221534720", "1717410430709663616"):
            continue
        for model in (
            "bj_photogeo",
            "spec_plx",
            "spdist_only",
            "dwarf_only",
            "gspphot",
        ):
            sample = draws(gs[s["source_id"]], s, model)
            if sample is None:
                continue
            _, x, v = sample
            phi, flag = eval_potential_c(
                pot, np.hypot(x[0], x[1]) / pot._ro, x[2] / pot._ro
            )
            assert flag == 0
            ves = np.sqrt(-2 * phi) * pot._vo
            vel = np.linalg.norm(v, axis=0)
            out.append(
                dict(
                    source_id=s["source_id"],
                    model=model,
                    n=len(vel),
                    vesc_16=np.quantile(ves, 0.16),
                    vesc_50=np.median(ves),
                    vesc_84=np.quantile(ves, 0.84),
                    p_unbound=np.mean(vel > ves),
                )
            )
    save("heavy_halo.csv", out)
