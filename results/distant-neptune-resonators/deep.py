"""Measure the duration and clone dependence of the VZ98 23:2 episode."""

import csv
from concurrent.futures import ProcessPoolExecutor
from dynamics import WORK
from clones import run

r = {"designation": "445473 (2010 VZ98)", "spkid": "20445473", "p": 23, "q": 2}
if __name__ == "__main__":
    tasks = [
        (64, 0.5, 1000000, 100),
        (64, -0.5, -1000000, 100),
        (8, 0.5, 10000000, 200),
        (8, -0.5, -10000000, 200),
        (8, 0.25, 1000000, 100),
        (8, -0.25, -1000000, 100),
    ]
    with ProcessPoolExecutor(6) as pool:
        for config, batch in zip(tasks, pool.map(run, [(r, *c) for c in tasks])):
            path = WORK / ("deep_%d_%s_%d.csv" % (config[0], config[1], config[2]))
            with path.open("w") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        "designation",
                        "spkid",
                        "p",
                        "q",
                        "clone",
                        "years",
                        "step_years",
                        "cadence_years",
                        "arc_deg",
                        "center_deg",
                        "unwrapped_range_deg",
                    ]
                )
                writer.writerows(batch)
            print(
                path.name, sum(row[-1] < 355 for row in batch), len(batch), flush=True
            )
