"""Check the saved public origin catalogs for names and comparable heliocentric orbits."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

SOURCES = {
    "tab1.tsv": "https://planet4589.org/space/deepcat/tab1.tsv",
    "tab2.tsv": "https://planet4589.org/space/deepcat/tab2.tsv",
    "sat_id2.htm": "https://projectpluto.com/sat_id2.htm",
}
ALIASES = ["2007 VE138", "2026 VV15", "7V30455", "LVg0046"]


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("directory",type=Path,help="directory containing the three source snapshots")
    args = cli.parse_args()
    result = {"retrieved_catalogs": {}, "search_aliases": ALIASES}
    for name, url in SOURCES.items():
        raw = (args.directory/name).read_bytes()
        text = raw.decode(errors="replace").casefold()
        result["retrieved_catalogs"][name] = {
            "url": url, "sha256": hashlib.sha256(raw).hexdigest(),
            "exact_string_counts": {alias: text.count(alias.casefold()) for alias in ALIASES},
        }
    with (args.directory/"tab2.tsv").open() as stream:
        phases = [row for row in csv.reader(stream, delimiter="\t")
                  if len(row)>=9 and row[3].strip()=="Sun" and
                  re.match(r"\d{4}",row[4].strip()) and int(row[4].strip()[:4])<=2007]
    numeric = []
    for row in phases:
        values = re.fullmatch(r"\s*([\d.]+)\??\s*x\s*([\d.]+)\??\s*AU\s*x\s*([\d.]+)\??\s*", row[-1])
        if values:
            q, Q, inc = map(float,values.groups())
            numeric.append((row[1].strip(),q,Q,inc))
    result["sun_centered_rows_with_start_year_through_2007"] = len(phases)
    result["rows_with_numeric_q_Q_i_and_bound_orbit"] = len(numeric)
    result["comparison_box"] = {"q_au": [0.69,0.89], "Q_au": [1.73,2.13],
                                "inclination_deg": [14.3,24.3]}
    result["matching_phases_in_box"] = [row for row in numeric if .69<=row[1]<=.89
                                        and 1.73<=row[2]<=2.13 and 14.3<=row[3]<=24.3]
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
