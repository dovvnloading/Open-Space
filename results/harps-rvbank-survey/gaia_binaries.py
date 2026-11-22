"""Read independently fitted DR3 orbits for RV survivors with binary flags."""
import argparse
import re
from pathlib import Path

import pandas as pd
from runtime import require_survey_workflow
from vet import GAIA, Sources


def run(evidence, out, cache):
    require_survey_workflow()
    vet = pd.read_csv(evidence / "vetting.csv", dtype=str, keep_default_na=False)
    selected = vet[vet.rv_status.eq("pass") & ~vet.binary_check.isin(["no catalog flag", "incomplete"])]
    hosts = selected[["star", "gaia_dr3_id"]].drop_duplicates()
    ids = sorted(i for i in set(hosts.gaia_dr3_id) if re.fullmatch(r"[0-9]{1,19}", i))
    db = Sources(cache, out)
    for table, columns in {
        "nss_two_body_orbit": ["source_id", "nss_solution_type", "period", "period_error",
                               "eccentricity", "semi_amplitude_primary", "semi_amplitude_primary_error"],
        "binary_masses": ["source_id", "m1", "m1_lower", "m1_upper", "m2", "m2_lower", "m2_upper"],
    }.items():
        sql = f"SELECT {','.join(columns)} FROM gaiadr3.{table} WHERE source_id IN ({','.join(ids)})"
        rows = db.get("Gaia DR3 " + table, GAIA, sql) if ids else None
        if rows is None or not set(columns).issubset(rows):
            rows = pd.DataFrame(columns=columns)
            available = False
        else:
            available = True
        joined = hosts.merge(rows, left_on="gaia_dr3_id", right_on="source_id", how="left", validate="many_to_many")
        joined["query_outcome"] = "no published row for ID" if available else "unavailable"
        joined.loc[joined.source_id.notna(), "query_outcome"] = "published DR3 row"
        joined.loc[joined.gaia_dr3_id.eq(""), "query_outcome"] = "no exact Gaia DR3 identifier"
        joined.to_csv(out / f"gaia-{table}.csv", index=False)
    (out / "catalog-inputs.csv").rename(out / "gaia-inputs.csv")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["evidence", "out", "cache"]:
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    run(a.evidence, a.out, a.cache)
