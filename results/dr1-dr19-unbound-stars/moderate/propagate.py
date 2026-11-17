import glob, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from infer import summarize
from make_queue import rows, save

if __name__ == "__main__":
    gs = {r["source_id"]: r for p in glob.glob("gaia/rows_*.csv") for r in rows(p)}
    out = []
    for s in rows("queue.csv"):
        out += summarize(s, gs[s["source_id"]])
    save("inference.csv", out)
