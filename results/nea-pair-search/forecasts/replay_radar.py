"""Compute reception-time, monostatic DSS-14 geometric delay and Doppler model differences."""
import argparse
from datetime import datetime,timedelta
from importlib import util
import json
import os
from pathlib import Path
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
LIGHT_KM_S = 299792.458
FREQ_HZ = 8.56e9
TIMES = ["2026-11-15T16:54:00Z", "2026-11-15T20:02:30Z", "2026-11-15T23:11:00Z"]


def uncertainty(engine, state, cov):
    errors = np.sqrt(np.diag(cov))
    weights, axes = np.linalg.eigh(cov / np.outer(errors, errors))
    differences = []
    for weight, axis in zip(weights, axes.T):
        d = np.sqrt(max(0, weight)) * errors * axis
        differences.append((engine(state + d) - engine(state - d)) / 2)
    return np.sqrt(np.sum(np.square(differences), axis=0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("radar_dir", type=Path)
    parser.add_argument("--orbits", type=Path, default=HERE.parent / "orbits")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    build = args.build.resolve()
    orbit_dir = args.orbits.resolve()
    spec = util.spec_from_file_location("orbit_replay", orbit_dir / "replay.py")
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.HERE = HERE
    os.environ["FO_ACCEL"] = "1"
    gravity = json.loads((orbit_dir / "fits/VE138.json").read_text())
    acceleration = json.loads((HERE / "force_states.json").read_text())
    gs = np.array(gravity["full"]["state"])
    states = {k: np.array(v) for k,v in acceleration["model_states"].items()}
    rows = []
    for utc in TIMES:
        t = datetime.fromisoformat(utc.replace("Z", "+00:00"))
        for offset in (-60, 0, 60):
            u = (t + timedelta(seconds=offset)).strftime("%Y-%m-%dT%H:%M:%SZ")
            rows.append(dict(utc=u, site="253", ra_deg="282.8", dec_deg="42.3",
                             catalog="Gaia3", mode="CCD"))
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder)
        wrapper = path / "find_orb"
        wrapper.mkdir()
        source = build / "find_orb"
        for f in source.iterdir():
            if f.is_file() and f.name != "fo":
                (wrapper / f.name).symlink_to(f)
        (wrapper / "fo").symlink_to(args.radar_dir.resolve() / "fo")
        (path / "de440.bin").symlink_to(build / "de440.bin")
        work = path / "run"
        work.mkdir()
        engine = module.Forward(path, work, rows, "2026 VV15", str(gravity["epoch_jd_tdb"]), gs)
        try:
            gr = engine(np.r_[gs,0.,0.])
            acc = {k:engine(x) for k,x in states.items()}
            cs = uncertainty(engine, states["A12"], np.asarray(acceleration["A12_site_utc_date_cluster_covariance"]))
            gcov = np.asarray(gravity["full"]["covariance"]["site_utc_date_clustered"])
            sg = uncertainty(lambda x: engine(np.r_[x, 0., 0.]), gs, gcov)
        finally:
            engine.close()
    res = []
    for i,utc in enumerate(TIMES):
        j = i*3 + 1
        models = {}
        for k, v in acc.items():
            change_delay, change_freq = v[j] - gr[j]
            delta_range = change_delay * LIGHT_KM_S / 2
            delta_speed = -change_freq * (LIGHT_KM_S * 1000) / (2*FREQ_HZ)
            models[k] = dict(delta_two_way_delay_us=round(change_delay*1e6,3),
                             delta_doppler_hz=round(change_freq,4),
                             delta_delay_equivalent_range_km=round(delta_range,5),
                             delta_doppler_equivalent_range_rate_ms=round(delta_speed,5))
        # Independent sign/unit check: radar Doppler is -f * d(round-trip delay)/dt.
        delta_delays = acc["A12"][j-1:j+2,0] - gr[j-1:j+2,0]
        fd = -FREQ_HZ * (delta_delays[-1] - delta_delays[0]) / 120
        if abs(fd - models["A12"]["delta_doppler_hz"]) > 0.4:
            raise ValueError("Doppler finite difference disagrees with two-leg projection")
        res.append(dict(utc_receive=utc, gravity_two_way_delay_s=round(gr[j,0],9),
                        gravity_echo_minus_transmit_hz=round(gr[j,1],5),
                        gravity_cluster_sigma_delay_us_doppler_hz=[round(sg[j,0]*1e6,2),round(sg[j,1],3)],
                        A12_cluster_sigma_delay_us_doppler_hz=[round(cs[j,0]*1e6,2),round(cs[j,1],3)],
                        A12_finite_difference_delta_doppler_hz=round(fd,4),
                        differences_from_gravity=models))
    output = {"site":"253, Goldstone DSS-14 monostatic, 8560 MHz",
              "status":"Hypothetical geometry. DSS-14 operations were suspended for modernization. These are calculated center-of-mass echoes with no measured radar return.",
              "utc_convention":"Reception epoch. Find_Orb light-time-corrects the target on receive, iterates the emitting station's retarded location, adds the two light times and computes the two-leg relativistic line-of-sight frequency factor. Doppler is modeled receive minus transmit frequency, and differences are force minus gravity. The response omits shape, rotation, media, instrumental and station-frequency biases.",
              "epochs":res}
    if args.output:
        args.output.write_text(json.dumps(output,indent=2)+"\n")
    else:
        previous = json.loads((HERE / "radar.json").read_text())
        for old,new in zip(previous["epochs"],res,strict=True):
            for k,v in new["differences_from_gravity"].items():
                ref = old["differences_from_gravity"][k]
                if abs(v["delta_two_way_delay_us"]-ref["delta_two_way_delay_us"])>0.1 or abs(v["delta_doppler_hz"]-ref["delta_doppler_hz"])>0.02:
                    raise RuntimeError("Replay mismatch at "+new["utc_receive"])
            print(new["utc_receive"],new["differences_from_gravity"]["A12"])


if __name__ == "__main__":
    main()
