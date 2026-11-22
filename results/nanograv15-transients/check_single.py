"""Small deterministic linear-algebra and injection checks for the manual workflow."""
import numpy as np

from search import (Projection, delete_visit, epoch_grid, measurable, score,
                    sensitivity, waveform, wideband_dm_check)


def example():
    rng = np.random.default_rng(73258)
    nights = 56200 + np.arange(96)*8 + rng.uniform(.1, .3, 96)
    t = np.repeat(nights, 2)
    f = np.tile([860., 1460.], len(nights))
    backend = np.where(t > np.median(t), "system-b", "system-a")
    groups = np.repeat(np.arange(len(nights)), 2)
    err = rng.uniform(.6, 1.3, len(t))
    jitter = .7 * np.ones(len(t))
    noise = rng.normal(size=len(t))*err + rng.normal(size=len(nights))[groups]*jitter
    epoch = float(np.mean(t))
    h = np.exp(-.5*((t-epoch)/14.)**2)
    data = {"t": t, "freq": f, "y": 120*h + noise, "backend": backend, "err": err,
            "design": np.sin((t-t[0])/84)[:, None], "ra": .4, "dec": -.8,
            "ecorr_group": groups, "ecorr_us": jitter, "dm_retained": True,
            "dm": .001 * h + rng.normal(size=len(t))*1.e-6,
            "dm_err": np.full(len(t), 1.e-6)}
    return data, epoch


def main():
    data, epoch = example()
    p = Projection(data, fit_soft=False)
    n = len(p.t)
    same_block = data["ecorr_group"][:, None] == data["ecorr_group"]
    covariance = np.diag(p.err**2) + same_block*np.outer(data["ecorr_us"], data["ecorr_us"])
    w = p.weight(np.eye(n))
    np.testing.assert_allclose(w @ covariance @ w.T, np.eye(n), atol=1.e-12)
    np.testing.assert_allclose(p.q.T @ p.q, np.eye(p.rank), atol=1.e-11)
    np.testing.assert_allclose(p.whiten(data["design"]), 0., atol=1.e-10)
    phase = 2*np.pi*(p.t - p.t.mean())/np.ptp(p.t)
    np.testing.assert_allclose(p.whiten(p.chrom*np.cos(phase)), 0., atol=1.e-10)
    raw = waveform(p, "achromatic_pulse", epoch, 14.).ravel()
    _, _, _, amp, _ = score(p, p.whiten(raw[:, None]))
    assert abs(amp[0]-120.) < 5., amp

    # Compare the efficient deletion identity to an explicit joint full-model refit.
    check = delete_visit(p, raw)
    assert check["status"] == "tested", check
    keep = np.floor(p.t) != check["mjd_day"]
    full = np.column_stack((p.q[keep], p.weight(raw)[keep]))
    w_y = p.weight(p.y)[keep]
    coeff, _, rank, _ = np.linalg.lstsq(full, w_y, rcond=None)
    assert rank == full.shape[1]
    scale = np.linalg.norm(w_y - full @ coeff)/np.sqrt(keep.sum()-rank)
    se = np.sqrt(np.linalg.inv(full.T @ full)[-1, -1])*scale
    np.testing.assert_allclose([check["amplitude_us"], check["t"]],
                               [coeff[-1], coeff[-1]/se], rtol=1.e-9)

    dm = wideband_dm_check(p, {"template": "dm_pulse", "epoch_mjd": epoch,
                               "width_d": 14., "amplitude_us": 2.})
    assert dm["status"] == "tested" and dm["directional_p"] < .001, dm
    np.testing.assert_allclose(dm["delta_dm_pc_cm3"], .001, rtol=.01)

    grid = epoch_grid(p.t)
    bank = waveform(p, "achromatic_pulse", grid, 14.)
    norm, unit, _, _, _ = score(p, p.whiten(bank))
    usable = measurable(p, bank, norm)
    result = sensitivity(p, "SYNTHETIC", "achromatic_pulse", 14., grid, norm, unit, usable, 4.)
    level = result["upper_limit_us_if_no_primary_crossing"]
    assert level is not None, result
    record = next(r for r in result["efficiency_nearest_matching_template"] if r["amplitude_us"] == level)
    assert record["recovered"] >= result["upper_order_rank"], result
    print("Passed release-covariance, nuisance projection, fixed-amplitude, block deletion, direct DM and injection checks.")


if __name__ == "__main__":
    main()
