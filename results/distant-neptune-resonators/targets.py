"""Propagate independent orbit fits for the two strongest previously unlisted objects."""

from concurrent.futures import ProcessPoolExecutor
import csv
from dynamics import WORK
from clones import run

if __name__ == "__main__":
    targets = [("2017 DO171", "54464562", 8, 1), ("2022 GV6", "54276384", 29, 3)]
    settings = [
        (96, 0.5, 400000, 100),
        (96, -0.5, -400000, 100),
        (48, 0.5, 2000000, 200),
        (48, -0.5, -2000000, 200),
        (8, 0.25, 400000, 100),
        (8, -0.25, -400000, 100),
    ]
    tasks = []
    for name, id_, p, q in targets:
        record = dict(designation=name, spkid=id_, p=p, q=q)
        tasks.extend((record, *s) for s in settings)
    with ProcessPoolExecutor(6) as pool:
        for task, batch in zip(tasks, pool.map(run, tasks)):
            r, n, step, years, interval = task
            path = WORK / ("target_%s_%d_%d.csv" % (r["spkid"], n, years))
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
            print(path.name, sum(row[-1] < 355 for row in batch[1:]), n, flush=True)
