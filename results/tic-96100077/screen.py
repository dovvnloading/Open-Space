"""Recreate both the known-system control count and the ranked combined-orbit list."""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
population = pd.read_csv(ROOT / "data" / "gaia.csv", dtype={"source_id": str})
table = (ROOT / "data" / "known.tsv").read_text().splitlines()
rows = [line.split("\t") for line in table if line.split("\t")[0].strip().isdigit()]
known = {row[2].strip() for row in rows}
if len(known) != 403:
    raise ValueError("Incomplete reference catalog")
population["in_czavalinga"] = population.source_id.isin(known)
population["photometric_period_days"] = 1 / population.frequency
follow = population[
    population.nss_solution_type.eq("AstroSpectroSB1") & ~population.in_czavalinga
]
follow = follow.sort_values("global_ranking", ascending=False)
follow.to_csv(ROOT / "outputs" / "ranked-combined-solutions.csv", index=False)
print(population.groupby(["in_czavalinga", "nss_solution_type"]).size().to_string())
