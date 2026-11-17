"""Refresh catalogue queries for Ross 154."""
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json

# Catalogue query machinery is shared with the adjacent MCC124 literature folder.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stellar"))
import mcc124_queries as queries
queries.OUT = Path(__file__).resolve().parent / "catalogue"
queries.OUT.mkdir(exist_ok=True)
gaia = queries.GAIA.replace("153.57985,21.07495,0.00834", "282.4606,-23.8376,0.014")
exo = """SELECT DISTINCT hostname,pl_name FROM ps WHERE hip_name='HIP 92403' OR gaia_dr3_id='Gaia DR3 4075141768785646848' OR tic_id='TIC 224752354' OR hostname IN ('Ross 154','GJ 729','HIP 92403','V1216 Sgr')"""
with ThreadPoolExecutor(max_workers=2) as pool:
    data = dict(pool.map(queries.get_one, [
        ("gaia_dr3", "https://gea.esac.esa.int/tap-server/tap/sync", gaia),
        ("confirmed_planets", "https://exoplanetarchive.ipac.caltech.edu/TAP/sync", exo)]))
(queries.OUT / "query.json").write_text(json.dumps(data, indent=2) + "\n")
