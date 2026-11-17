"""Build the targeted SDSS comparison set, including a known HVS and RRc control."""

import numpy as np
from zero_point import zpt
from make_queue import rows, save, number

# Repeat-visit leads from the DR19 screen, one literature positive and one pulsation control.
IDS = [
    "4406622293782568192",
    "1629431706439574528",
    "869830353178482560",
    "2641617578115676160",
    "4337459232822884992",
    "1170190610327078144",
    "4433836684037871232",
]


def run():
    zpt.load_tables()
    gs = {r["source_id"]: r for r in rows("sdss/gaia_tail.csv")}
    out = []
    for sid in IDS:
        hit = next(
            r for r in rows("sdss/prefilter_tail.csv") if r["gaia_dr3_source_id"] == sid
        )
        g = gs[sid]
        zp = zpt.get_zpt(
            *[
                np.array([number(g[k])])
                for k in (
                    "phot_g_mean_mag",
                    "nu_eff_used_in_astrometry",
                    "pseudocolour",
                    "ecl_lat",
                    "astrometric_params_solved",
                )
            ],
            _warnings=False,
        )[0]
        s = dict(
            source_id=sid,
            targetid=hit["sdss_id"],
            channel="sdss",
            rv=hit["v_rad"],
            e_rv=hit["e_v_rad"],
            zp=zp,
            **{
                k: v
                for k, v in hit.items()
                if k.startswith(("r_med_", "r_lo_", "r_hi_"))
            },
        )
        if not np.isfinite(zp):
            args = [
                np.array([number(g[k])])
                for k in (
                    "phot_g_mean_mag",
                    "nu_eff_used_in_astrometry",
                    "pseudocolour",
                    "ecl_lat",
                    "astrometric_params_solved",
                )
            ]
            s["zp"] = zpt.get_zpt(*args, _warnings=True)[0]
            s["zp_flag"] = "colour outside calibration"
        if sid in {"1629431706439574528", "869830353178482560", "2641617578115676160"}:
            s["channel"] = "blend_control"
        if sid == "4406622293782568192":
            s.update(spec_d=2.733, spec_e=2.733 * np.sqrt(np.expm1(0.420**2)))
        if sid == "1170190610327078144":
            s.update(
                spec_d=4.431,
                spec_e=0.48,
                rv=67.5897,
                e_rv=2.3706,
                channel="RRc_control",
            )
        if sid == "4337459232822884992":
            s["channel"] = "known_control"
        if number(s["e_rv"]) <= 0:
            good = [
                number(r.get("e_v_rel"))
                for r in rows("sdss/mwm_apogee_allvisit_tail.csv")
                if r["sdss_id"] == hit["sdss_id"] and number(r.get("e_v_rel")) > 0
            ]
            s["e_rv"] = (
                max(1.0, float(np.median(good)))
                if good
                else max(1.0, number(hit["std_v_rad"]))
            )
        out.append(s)
    save("sdss_queue.csv", out)


if __name__ == "__main__":
    run()
