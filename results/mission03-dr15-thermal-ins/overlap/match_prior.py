"""Audit source-level CSVs against published candidate and repository positions."""

import argparse
import csv
import math
from pathlib import Path


def separation(ra1, dec1, ra2, dec2):
    """Great-circle distance in arcsec, including wrap at RA = 0."""
    dra = math.radians(ra2 - ra1)
    d1, d2 = math.radians(dec1), math.radians(dec2)
    h = (
        math.sin((d2 - d1) / 2) ** 2
        + math.cos(d1) * math.cos(d2) * math.sin(dra / 2) ** 2
    )
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, h)))) * 3600


def run(inputs, priors, output):
    with priors.open(newline="") as stream:
        reference = list(csv.DictReader(stream))
    fields = [
        "SRCID",
        "IAUNAME",
        "review_count",
        "prior",
        "separation_arcsec",
        "review_radius_arcsec",
        "within_review_radius",
        "position_source",
    ]
    with inputs.open(newline="") as src, output.open("w", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for raw in csv.DictReader(src):
            row = {key.upper(): value for key, value in raw.items()}
            ra, dec = float(row["RA"]), float(row["DEC"])
            nearest = []
            for prior in reference:
                sep = separation(
                    ra, dec, float(prior["ra_deg"]), float(prior["dec_deg"])
                )
                # Triage envelope for aliases, catalogue revisions and older positions.
                # Final identity requires checking each counterpart's actual astrometry.
                radius = 90.0 if "RXS J" in prior["designation"] else 30.0
                nearest.append((sep, radius, prior))
            nearest.sort(key=lambda item: item[0])
            matches = [item for item in nearest if item[0] <= item[1]]
            for sep, radius, prior in matches or nearest[:1]:
                writer.writerow(
                    dict(
                        SRCID=row["SRCID"],
                        IAUNAME=row.get("IAUNAME") or row["NAME"],
                        review_count=len(matches),
                        prior=prior["designation"],
                        separation_arcsec=f"{sep:.3f}",
                        review_radius_arcsec=radius,
                        within_review_radius=int(sep <= radius),
                        position_source=prior["position_source"],
                    )
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_csv", type=Path)
    parser.add_argument("prior_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    run(args.source_csv, args.prior_csv, args.output_csv)
