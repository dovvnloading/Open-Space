"""Run the complementary short/missing-reported-arc stratum with the shared fitter."""
import importlib.util
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("ground_census",ROOT.parent/"census/replay.py")
census=importlib.util.module_from_spec(spec)
# ProcessPoolExecutor pickles worker functions by module name in forked children.
sys.modules[spec.name]=census
spec.loader.exec_module(census)

if __name__=="__main__":
    census.CACHE=ROOT/"cache"
    census.OUT=ROOT/"evidence"
    census.CACHE.mkdir(exist_ok=True)
    census.OUT.mkdir(exist_ok=True)
    census.acquire.CACHE=census.fit.CACHE=census.fit.DATA=census.CACHE
    census.main(Path(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),short=True)
