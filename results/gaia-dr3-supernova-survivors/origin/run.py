"""Escape speeds and Milky Way disk crossings for spectroscopic controls."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CONTROL = HERE.parent / "control"
PRIOR = HERE.parents[1] / "dr1-dr19-unbound-stars"
sys.path.insert(0, str(PRIOR))
from fetch_gaia import fetch
from infer import draws, escape
from make_queue import rows, save
from orbits import integrate

FOCUS = [["D6-1", "D6-2", "D6-3", "J1109+0001"], ["LP 40-365", "J1240+6710", "US 708"]]


def indexed(path):
    return {r["source_id"]: r for r in rows(path)}


def main():
    controls, gaia, receipts = [], {}, []
    for batch, names in enumerate(FOCUS):
        group = [r for r in rows(CONTROL/"controls.csv") if r["name"] in names]
        ids = sorted([r["source_id"] for r in group], key=int)
        receipt = fetch((batch, ids, HERE))
        data = indexed(HERE/receipt["file"])
        if set(data) != set(ids):
            raise ValueError("Gaia orbit input is incomplete")
        controls.extend(group)
        gaia.update(data)
        receipts.append(receipt)
    (HERE/"gaia_receipt.json").write_text(json.dumps(receipts, indent=2)+"\n")
    kin = indexed(CONTROL/"published_kinematics.csv")
    bj = indexed(CONTROL/"bj.csv")
    calibration = indexed(CONTROL/"recovery.csv")
    states, crossings = [], []
    for star in controls:
        sid = star["source_id"]
        radial = kin[sid]
        s = dict(bj[sid], rv=radial["rv_kms"], e_rv=radial["rv_err_kms"],
                 rv_floor=1.0, zp=calibration[sid]["zp_mas"])
        models = ["bj_geo", "bj_photogeo", "edsd_plx"]
        if star["name"] == "US 708":
            # Geier et al. spectroscopic distance used in El-Badry et al. 2026.
            s.update(spec_d=8.5, spec_e=1.0)
            models.append("spec_plx")
        for model in models:
            sample = draws(gaia[sid], s, model)
            if sample is None:
                continue
            d, x, v = sample
            speed = np.linalg.norm(v, axis=0)
            mc, gala = escape(x)
            result = dict(source_id=sid, name=star["name"], distance_model=model, draws=len(d),
                          rv_kms=s["rv"], e_rv_kms=s["e_rv"],
                          p_escape_mc17=float(np.mean(speed>mc)),
                          p_escape_gala=float(np.mean(speed>gala)))
            for key, data in [("distance_kpc",d), ("vgal_kms",speed), ("vesc_mc17_kms",mc),
                              ("vesc_gala_kms",gala)]:
                for pct, vq in zip((16,50,84), np.percentile(data,[16,50,84])):
                    result[f"{key}_{pct}"] = round(float(vq), 3)
            states.append(result)
            for dt in (0.25,0.125):
                cross = integrate(s, gaia[sid], model, dt=dt, n=256, duration=50.0)
                crossings.append(dict(name=star["name"], **cross))
    save(HERE/"states.csv",states)
    save(HERE/"crossings.csv",crossings)


if __name__ == "__main__":
    main()
