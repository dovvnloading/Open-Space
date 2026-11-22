"""Bridge the published TEMPO2 EPTA DR2 timing solution to PINT 1.1.8."""
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path
from tempfile import TemporaryDirectory

import astropy.units as u
import numpy as np
from pint.fitter import Fitter
from pint.models import get_model_and_toas
from pint.observatory import get_observatory
from pint.residuals import Residuals
from pint.toa import _parse_TOA_line, _toa_format, get_TOAs_array


TIM_COMMANDS = set("INCLUDE INFO FORMAT EFAC EQUAD EMIN EMAX FMIN FMAX MODE PHASE TIME JUMP".split())
UNSUPPORTED_COMMANDS = set("EFLOOR ESET SIGMA GLOBAL_EFAC T2EFAC T2EQUAD".split())


def toa_columns(fields):
    """The FORMAT 1 reader needs filename, frequency, site MJD, error and site."""
    if len(fields) < 5:
        return False
    try:
        float(fields[1])
        float(fields[3])
    except ValueError:
        return False
    return True


def is_flag(token):
    return len(token) >= 2 and token[0] == "-" and (token[1].isalpha() or token[1] == "_")


def normalize_flags(tokens, counts, model_tokens):
    pairs = []
    offset = Decimal(0)
    i = 0
    while i < len(tokens):
        k = tokens[i]
        i += 1
        if not is_flag(k):
            # TEMPO2 scans for flag identifiers; bare trailing notes are not flags.
            counts["nonflag_tokens_ignored"] += 1
            continue
        if i == len(tokens) or is_flag(tokens[i]):
            if k.lower() in ("-addsat", "-to", "-pn", "-phase") or k in model_tokens:
                raise ValueError(f"missing value for timing flag or published model selector {k}")
            counts["valueless_flags_ignored"] += 1
            continue
        v = tokens[i]
        i += 1
        label = k.lstrip("-")
        if label in ("error", "freq", "scale", "MJD", "flags", "obs", "name"):
            if k in model_tokens:
                raise ValueError(f"PINT-reserved flag used in timing model: {k}")
            pairs.extend(("-tempo2_meta_"+label, v))
            counts["renamed_pint_reserved_metadata_flags"] += 1
        else:
            pairs.extend((k, v))
        if k == "-addsat":
            # TEMPO2 adds seconds to the *site* time before clock corrections.
            # PINT otherwise treats this flag as inert metadata.
            offset += Decimal(v)
            counts["addsat_flags_applied"] += 1
    return pairs, offset


def stage_tim(src, work, record, model_tokens):
    lines = src.read_text().splitlines(keepends=True)
    parsed = [s.split() for s in lines]
    declared = any(f[:2] == ["FORMAT","1"] for f in parsed)
    modern = declared or any(toa_columns(f) and "." in f[2] and len(f[2].split(".")[0]) == 5
                             for f in parsed)
    counts = Counter()
    file_id = str(src.relative_to(work))
    record["tim_file_counts"][file_id] = counts
    transformed = []
    if modern and not declared:
        transformed.append("FORMAT 1\n")
        counts["inferred_format1_header"] += 1
    skip = False
    profile_dir = None
    for line_no, (line, fields) in enumerate(zip(lines, parsed), 1):
        if not fields or line.startswith("#"):
            transformed.append(line)
            continue
        word = fields[0].upper()
        if (modern and line.startswith("C")) or (not modern and fields[0] in ("C","c")):
            record["normalized_comments"] += 1
            counts["tempo2_comments"] += 1
            continue
        if word == "NOSKIP":
            skip = False
            continue
        if skip:
            counts["skipped_lines"] += 1
            continue
        if word == "SKIP":
            skip = True
            continue
        if word == "END":
            counts["lines_after_local_end"] = len(lines)-line_no
            break
        if word == "PROFILE_DIR":
            if len(fields) < 2:
                raise ValueError(f"{file_id}:{line_no}: PROFILE_DIR has no value")
            profile_dir = fields[1]
            counts["profile_directory_directives"] += 1
            continue
        deleted = modern and line.startswith("I") and toa_columns(fields[1:])
        if deleted:
            # INCLUDE/INFO commands do not have the additional filename and
            # numeric TOA columns of a TEMPO2 I-prefixed deleted observation.
            counts["tempo2_I_deleted_toas"] += 1
            continue
        if word in TIM_COMMANDS and not toa_columns(fields):
            transformed.append(line)
            continue
        if word in UNSUPPORTED_COMMANDS:
            raise ValueError(f"{file_id}:{line_no}: unsupported timing-error directive {word}")
        if not modern:
            transformed.append(line)
            continue
        if not toa_columns(fields) or line.startswith("I"):
            # TEMPO2's sscanf needs five parsed fields (six for a deleted row).
            # Invalid nondirective lines are ignored, including free-form notes.
            counts["incomplete_nondirective_lines"] += 1
            record["ignored_line_locations"].setdefault(file_id, []).append(line_no)
            continue
        try:
            mjd_value = Decimal(fields[2])
            pairs, offset = normalize_flags(fields[5:], counts, model_tokens)
            if offset:
                with localcontext() as context:
                    context.prec = 45
                    mjd_value += offset/Decimal(86400)
            if offset or "e" in fields[2].lower():
                fields[2] = format(mjd_value,"f")
            mjd = float(mjd_value)
            if not np.isfinite(mjd) or not np.isfinite(float(fields[1])) or not np.isfinite(float(fields[3])):
                raise ValueError("nonfinite timing measurement")
        except (ArithmeticError, ValueError) as e:
            raise ValueError(f"{file_id}:{line_no}: invalid listed TOA: {e}") from e
        site = fields[4]
        sites = record["clock_valid_mjd_by_site"]
        if site not in sites:
            sites[site] = site_clock_ranges(site)
        missing = [clock for clock, (start, end) in sites[site].items() if not start <= mjd <= end]
        if missing:
            identity = file_id + ":" + site + ":" + ",".join(missing)
            cut = record["clock_cuts"].setdefault(identity, {"n":0, "min_mjd":mjd, "max_mjd":mjd})
            cut["n"] += 1
            cut["min_mjd"] = min(mjd, cut["min_mjd"])
            cut["max_mjd"] = max(mjd, cut["max_mjd"])
            counts["outside_clock_range"] += 1
            continue
        # A tab after the filename prevents PINT from misclassifying a valid
        # one-character filename such as lowercase c as a Princeton record.
        name = fields[0] if profile_dir is None else profile_dir + "/" + fields[0]
        tail = "\t"+" ".join(fields[1:5]+pairs)+"\n"
        native = name+tail
        if _toa_format(native, "Tempo2") != "Tempo2":
            # PINT recognizes command *prefixes* even in valid TEMPO2 filenames.
            if "-name" in model_tokens or name in model_tokens:
                raise ValueError(f"{file_id}:{line_no}: ambiguous filename used by published timing model")
            native = "./"+name+tail
            counts["escaped_pint_command_prefix_filenames"] += 1
        try:
            _, parsed_toa = _parse_TOA_line(native, fmt="Tempo2")
            if parsed_toa["format"] != "Tempo2":
                raise ValueError("valid FORMAT 1 TOA misclassified")
        except (ValueError, IndexError, KeyError) as exc:
            raise ValueError(f"{file_id}:{line_no}: staged TOA incompatible with PINT: {exc}") from exc
        transformed.append(native)
        counts["retained_format1_toas"] += 1
    return "".join(transformed)


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
                                   "clock_valid_mjd_by_site": {}, "tim_file_counts": {},
                                   "ignored_line_locations": {}}
    sources = [tim] + [work / name for name in row["include_sha256"]]
    model_tokens = set(par.read_text().split())
    with TemporaryDirectory() as tmp:
        for src in sources:
            target = Path(tmp) / src.relative_to(work)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(stage_tim(src,work,record,model_tokens))
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
