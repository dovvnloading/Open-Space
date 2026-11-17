"""Recreate the additional screen and exclude sources already present in the first pass."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from screen_desi import screen
from make_queue import rows, save

if __name__ == "__main__":
    screen(sys.argv[1], 550)
    first = {r["SOURCE_ID"] for r in rows("../desi_screen.csv")}
    save(
        "desi_screen.csv",
        [r for r in rows("desi_screen.csv") if r["SOURCE_ID"] not in first],
    )
