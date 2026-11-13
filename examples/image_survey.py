"""Plan coordinate-only TESS image batches and classify their automatic follow-up needs."""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from open_space.tess_archive import atomic_bytes


def plan_grid(ra, dec, columns, rows, step_ra, step_dec, sectors, size, start, limit):
    if (
        not 0 <= ra < 360
        or not -90 < dec < 90
        or min(columns, rows, limit) < 1
        or start < 0
        or step_ra <= 0
        or step_dec <= 0
        or not 11 <= size <= 41
        or not 1 <= sectors[0] < sectors[1] <= 9999
    ):
        raise ValueError("Invalid field grid, sector pair, or page")
    if (
        start >= columns * rows
        or dec + (rows - 1) * step_dec >= 90
        or (columns - 1) * step_ra >= 360
    ):
        raise ValueError(
            "Grid exceeds declination limits, wraps more than once, or page is empty"
        )
    fields = []
    for n in range(start, min(start + limit, columns * rows)):
        y, x = divmod(n, columns)
        lon, lat = round((ra + x * step_ra) % 360, 9), round(dec + y * step_dec, 9)
        fields.append(
            {
                "field_id": f"grid-{lon:.9f}{lat:+.9f}-s{sectors[0]}-{sectors[1]}",
                "ra_deg": lon,
                "dec_deg": lat,
                "size_pixels": size,
                "role": "coordinate_grid_exploration",
                "products": [{"sector": s} for s in sectors],
            }
        )
    return {
        "schema_version": 1,
        "selection": {
            "frame": "ICRS",
            "ordering": "row-major starting at index zero",
            "origin_deg": [ra, dec],
            "shape": [rows, columns],
            "step_deg": [step_ra, step_dec],
            "start": start,
            "total_cells": columns * rows,
            "overlap": "Coordinate steps; neighboring cutouts can overlap.",
        },
        "fields": fields,
    }


def pin(plan, cache, workers):
    from open_space.ffi_archive import acquire

    if not 1 <= workers <= 16:
        raise ValueError("Choose 1 to 16 archive readers")

    def field_data(field):
        products = []
        for product in field["products"]:
            _, receipt = acquire(field, product, cache)
            products.append(
                {
                    "sector": product["sector"],
                    "science_sha256": receipt["science_sha256"],
                }
            )
        return {**field, "products": products}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        fields = list(pool.map(field_data, plan["fields"]))
    return {**plan, "fields": fields}


def triage(directory):
    """Select sources for another observation using only completed queue evidence."""
    registry = json.loads((directory / "registry.json").read_bytes())
    fields = []
    lookups = {}
    for job in registry["tasks"]:
        if job["state"] != "done" or job["kind"] not in {"ffi", "counterparts"}:
            continue
        data = json.loads((directory / job["evidence"]).read_bytes())
        if job["kind"] == "ffi":
            fields.append((data, job["id"]))
        elif "ra_deg" in job["payload"]:
            pos = job["payload"]
            lookups[(job["subject"], pos["ra_deg"], pos["dec_deg"])] = (data, job["id"])
    rows = []
    for field, task in fields:
        for source in field["sources"]:
            if not source["candidate"]:
                continue
            lookup = lookups.get(
                (field["field_id"], source["ra_deg"], source["dec_deg"])
            )
            data, catalog_task = lookup if lookup else ({}, None)
            catalog = data.get("variable_catalog", {})
            matches = []
            for match in catalog.get("vsx", {}).get("matches", []):
                p = match.get("period_days")
                if not p or match["separation_arcsec"] > 15:
                    continue
                if match.get("period_limit_flag") not in ("", None):
                    continue
                for scale in (0.5, 1, 2):
                    if all(
                        abs(part["period_days"] / (scale * p) - 1) <= 0.025
                        for part in [source["discovery"], source["validation"]]
                    ):
                        matches.append(
                            {
                                "counterpart": match["name"],
                                "period_days": p,
                                "factor": scale,
                                "separation_arcsec": match["separation_arcsec"],
                            }
                        )
            state = (
                "catalog_check_incomplete"
                if catalog.get("status") != "ok"
                else (
                    "catalog_period_association"
                    if matches
                    else "request_later_images_and_resolved_photometry"
                )
            )
            rows.append(
                {
                    "field": field["field_id"],
                    "source_id": source["source_id"],
                    "ra_deg": source["ra_deg"],
                    "dec_deg": source["dec_deg"],
                    "period_days": [
                        source[k]["period_days"] for k in ("discovery", "validation")
                    ],
                    "action": state,
                    "vsx_associations": matches,
                    "field_task": task,
                    "catalog_task": catalog_task,
                }
            )
    return {
        "method": "vsx-spatial-period-triage-v1",
        "radius_arcsec": 15,
        "fractional_tolerance": 0.025,
        "period_scales": [0.5, 1, 2],
        "rows": rows,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="command", required=True)
    g = subs.add_parser("plan")
    g.add_argument(
        "--origin", nargs=2, type=float, required=True, metavar=("RA", "DEC")
    )
    g.add_argument(
        "--shape", nargs=2, type=int, required=True, metavar=("COLUMNS", "ROWS")
    )
    g.add_argument(
        "--step", nargs=2, type=float, required=True, metavar=("RA_STEP", "DEC_STEP")
    )
    g.add_argument("--sectors", nargs=2, type=int, required=True)
    g.add_argument("--size-pixels", type=int, default=21)
    g.add_argument("--start", type=int, default=0)
    g.add_argument("--limit", type=int, default=32)
    f = subs.add_parser("pin")
    f.add_argument("plan", type=Path)
    f.add_argument("--cache", type=Path, required=True)
    f.add_argument("--workers", type=int, default=2)
    audit = subs.add_parser("triage")
    audit.add_argument("export_directory", type=Path)
    for cmd in (g, f, audit):
        cmd.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    if a.command == "plan":
        result = plan_grid(
            *a.origin, *a.shape, *a.step, a.sectors, a.size_pixels, a.start, a.limit
        )
    elif a.command == "pin":
        result = pin(json.loads(a.plan.read_bytes()), a.cache, a.workers)
    else:
        result = triage(a.export_directory)
    atomic_bytes(
        a.output, (json.dumps(result, indent=2, allow_nan=False) + "\n").encode()
    )


if __name__ == "__main__":
    main()
