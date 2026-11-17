from infer import load_inputs
from make_queue import save
from orbits import integrate

if __name__ == "__main__":
    stars, gs = load_inputs()
    s = next(r for r in stars if r["source_id"] == "1717410430709663616")
    g = gs[s["source_id"]]
    save(
        "convergence.csv",
        [
            integrate(s, g, model, dt)
            for model in ("spec_plx", "bj_photogeo")
            for dt in (0.25, 0.125)
        ],
    )
