"""Bridge the published TEMPO2 EPTA DR2 timing solution to PINT 1.1.8."""
from pathlib import Path
from tempfile import TemporaryDirectory

import astropy.units as u
import numpy as np
from pint.fitter import Fitter
from pint.models import get_model_and_toas
from pint.observatory import get_observatory
from pint.residuals import Residuals
from pint.toa import get_TOAs_array


def site_clock_ranges(code):
    """Use the exact clock chain and endpoint conventions from PINT's observatory definition."""
    site = get_observatory(code.upper())
    ranges = {}
    if hasattr(site, "_load_clock_corrections"):
        site._load_clock_corrections()
        for clock in site._clock:
            tabulated = getattr(clock, "clock_file", clock)
            if tabulated.valid_beyond_ends:
                continue
            if not len(tabulated.time):
                raise ValueError(f"no tabulated clock corrections for {code}: {clock.friendly_name}")
            name = Path(str(clock.friendly_name)).name
            ranges[name] = [float(tabulated.time[0].mjd), float(tabulated.time[-1].mjd)]
    return ranges


def load_epta(par, tim, work, row, eph, bipm):
    """Preserve TEMPO2 comments/INCLUDE paths and restrict TOAs to tabulated clocks."""
    # TEMPO2 readTimfile.C skips every physical line with C in column one.
    # PINT identifies only a subset of those lines, before its FORMAT 1 parser.
    # Reproduce the column-one rule in copies of the already-hashed inputs.
    record = row["epta_format"] = {"normalized_comments": 0, "clock_cuts": {},
                                   "clock_valid_mjd_by_site": {}}
    sources = [tim] + [work / name for name in row["include_sha256"]]
    sites = record["clock_valid_mjd_by_site"]
    with TemporaryDirectory() as tmp:
        for src in sources:
            target = Path(tmp) / src.relative_to(work)
            target.parent.mkdir(parents=True, exist_ok=True)
            transformed = []
            for line in src.read_text().splitlines(keepends=True):
                if line.startswith("C"):
                    record["normalized_comments"] += 1
                    transformed.append("#" + line[1:])
                    continue
                fields = line.split()
                if len(fields) >= 5 and not line.lstrip().startswith("#"):
                    try:
                        float(fields[1])
                        mjd = float(fields[2])
                        float(fields[3])
                    except ValueError:
                        pass  # Pass commands through to the original PINT parser.
                    else:
                        site = fields[4]
                        if site not in sites:
                            sites[site] = site_clock_ranges(site)
                        missing = [clock for clock, (start, end) in sites[site].items()
                                   if not start <= mjd <= end]
                        if missing:
                            key = str(src.relative_to(work)) + ":" + site + ":" + ",".join(missing)
                            cut = record["clock_cuts"].setdefault(key, {"n": 0, "min_mjd": mjd, "max_mjd": mjd})
                            cut["n"] += 1
                            cut["min_mjd"] = min(mjd, cut["min_mjd"])
                            cut["max_mjd"] = max(mjd, cut["max_mjd"])
                            transformed.append("# outside tabulated clock range; " + line)
                            continue
                transformed.append(line)
            target.write_text("".join(transformed))
        model, toas = get_model_and_toas(
            str(par), str(Path(tmp) / tim.relative_to(work)), usepickle=False,
            allow_name_mixing=True, allow_T2=True, allow_tcb=True, ephem=eph,
            include_bipm=True, bipm_version=bipm, limits="error",
            ell1h_shapiro="absorbed",
        )
    # PINT 1.1.8 otherwise uses the current default BIPM at the auxiliary
    # infinite-phase reference even when the real observations use BIPM2021.
    # Fill the same AbsPhase cache PINT uses, including its normal invalidation keys.
    ref = model.components["AbsPhase"]
    ref.tz_cache = get_TOAs_array(
        (ref.TZRMJD.quantity.jd1 - 2400000.5, ref.TZRMJD.quantity.jd2),
        obs=ref.TZRSITE.value, freqs=ref.TZRFRQ.quantity,
        include_bipm=True, bipm_version=bipm, ephem=eph,
        planets=toas.planets, tzr=True, limits="error",
    )
    ref.tz_hash = hash((ref.TZRMJD.value, ref.TZRSITE.value, ref.TZRFRQ.value))
    ref.tz_clkc_info = toas.clock_corr_info
    ref.tz_planets = toas.planets
    ref.tz_ephem = toas.ephem
    return model, toas


def refit_epta(model, toas, row):
    """Refit only released timing parameters at fixed pulse numbers after TCB conversion."""
    # The paper's TRES and NTOA must stay the comparison values. Fitter
    # updates both fields as part of its normal bookkeeping.
    original_tres = model.TRES.quantity
    original_ntoa = model.NTOA.value
    initial = Residuals(toas, model)
    initial_rms = float(initial.rms_weighted().to_value(u.us))
    record = row["epta_timing_refit"] = {
        "prefit_wrms_us": initial_rms, "initial_track_mode": initial.track_mode,
        "original_par_tres_us": None if original_tres is None else float(original_tres.to_value(u.us)),
        "original_par_ntoa": None if original_ntoa is None else int(original_ntoa),
        "initial_free_parameters": list(model.free_params),
    }
    if initial.track_mode != "use_pulse_numbers":
        # A cycle cannot be inferred from a badly smeared converted ephemeris.
        # Retain the phase connection from the supplied solution while refitting.
        if not np.isfinite(initial_rms) or initial_rms >= 0.05 * 1e6 / model.F0.value:
            raise ValueError("EPTA initial ephemeris is not sufficiently phase connected to refit")
        toas.compute_pulse_numbers(model)
        record["pulse_numbers"] = "nearest integer from the released converted model, fixed for all iterations"
    else:
        record["pulse_numbers"] = "published, fixed for all iterations"
    # A clock range cut can remove a backend completely. Its offset has zero
    # derivative on all retained TOAs and cannot be fit.
    missing = []
    for name in model.free_params:
        param = getattr(model, name)
        if hasattr(param, "select_toa_mask") and len(param.select_toa_mask(toas)) == 0:
            param.frozen = True
            missing.append(name)
    record["empty_masks_frozen"] = missing
    solver = Fitter.auto(toas, model, downhill=False, track_mode="use_pulse_numbers")
    solver.fit_toas(maxiter=3)
    result = solver.model
    # Keep explicit pulse-number tracking for all following residual calculations.
    if hasattr(result, "TRACK"):
        result.TRACK.value = "-2"
    post = Residuals(toas, result, track_mode="use_pulse_numbers")
    record.update(fitter=type(solver).__name__, iterations=3,
                  postfit_wrms_us=float(post.rms_weighted().to_value(u.us)),
                  max_abs_postfit_cycles=float(np.max(np.abs(post.phase_resids.value))))
    if not np.isfinite(record["max_abs_postfit_cycles"]) or record["max_abs_postfit_cycles"] >= .5:
        raise ValueError("EPTA refit leaves an ambiguous pulse count")
    return result
