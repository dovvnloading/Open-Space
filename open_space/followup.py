"""Execute the analyst's selected checks directly on the measured spectrum."""

import math
from statistics import median


def _usable(row):
    flux, weight = float(row["flux"]), float(row["ivar"])
    return (int(row["mask"]) == 0 and math.isfinite(flux)
            and math.isfinite(weight) and weight > 0)


def measure(requests, features, rows):
    results = []
    for request in requests:
        feature = features[int(request["feature_id"][1:]) - 1]
        center = feature["center_pixel"]
        core = rows[center - 1:center + 2]
        nearby = rows[max(0, center - 20):center + 21]
        result = {**request, "center_pixel": center, "observed_angstrom": feature["wavelength"]}
        if request["check"] == "mask_neighborhood":
            result.update(
                window_pixels=len(nearby), unusable_pixels=sum(not _usable(r) for r in nearby),
                core_mask=[int(r["mask"]) for r in core],
                core_ivar=[float(r["ivar"]) for r in core],
                core_or_mask=[int(r["or_mask"]) for r in core] if "or_mask" in rows[0] else None,
            )
        elif request["check"] == "continuum_sensitivity":
            # Two alternative sideband choices estimate the endpoint medians.
            # Interpolate to the center. Missing or masked support is inconclusive.
            result["fits"] = []
            for inner, outer in ((5, 12), (9, 20)):
                left = rows[max(0, center - outer):center - inner + 1]
                right = rows[center + inner:center + outer + 1]
                fit = {"inner_pixels": inner, "outer_pixels": outer}
                window = rows[max(0, center - outer):center + outer + 1]
                steps = [float(b["wavelength"]) - float(a["wavelength"]) for a, b in zip(window, window[1:])]
                if (len(window) != 2 * outer + 1 or not all(_usable(r) for r in window)
                        or max(steps) > 1.8 * median(steps)):
                    fit["status"] = "incomplete_coverage"
                else:
                    xl, xr = (median(float(r["wavelength"]) for r in side) for side in (left, right))
                    yl, yr = (median(float(r["flux"]) for r in side) for side in (left, right))
                    weight = sum(float(r["ivar"]) for r in core)
                    contrast = sum(float(r["ivar"]) * (float(r["flux"]) - yl -
                                   (yr - yl) * (float(r["wavelength"]) - xl) / (xr - xl)) for r in core)
                    fit.update(status="measured", core_contrast=contrast / weight)
                result["fits"].append(fit)
        else:
            raise ValueError("Unsupported follow-up measurement")
        results.append(result)
    return results
