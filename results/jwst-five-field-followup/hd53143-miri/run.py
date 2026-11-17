"""Run PR 142's F2550W search and injection experiment on each GO 8142 mosaic."""

import argparse
import sys
from pathlib import Path

from freeze_inputs import HERE, PREVIOUS, main as verify_inputs

sys.path.insert(0, str(PREVIOUS))
from reduce import analyze


def main(cache):
    rows = verify_inputs(cache, False)
    standards = [r for r in rows if r["role"] == "standard"]
    for science in (r for r in rows if r["role"] in ("roll", "combined")):
        visit = science["obs_id"].split("_")[0]
        analyze(
            cache,
            HERE / "measurements" / visit,
            [science, *standards],
            [("F2550W", 0.83)],
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    main(parser.parse_args().cache)
