"""Check MAST public and exclusive footprints, including other space-imaging bands."""

from urllib.request import Request, urlopen
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
import json
import pandas as pd


def query(spec):
    name, collection, instrument = spec
    filters = [{"paramName": "obs_collection", "values": [collection]}]
    if instrument:
        filters.append({"paramName": "instrument_name", "values": [instrument]})
    else:
        filters.append({"paramName": "t_min", "values": [{"min": 59580, "max": 61380}]})
    args = {
        "service": "Mast.Caom.Filtered.Position",
        "params": {
            "columns": "obs_id,t_min,t_max,filters,proposal_id,obs_collection,instrument_name,dataRights,calib_level",
            "filters": filters,
            "position": "150.122033, 2.323386208, 0.001",
        },
        "pagesize": 10000,
        "page": 1,
        "format": "json",
    }
    request = Request(
        "https://mast.stsci.edu/api/v0/invoke",
        data=urlencode({"request": json.dumps(args)}).encode(),
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "Mozilla/5.0 (JWST archive analysis)",
        },
    )
    result = json.load(urlopen(request, timeout=70))
    if result["status"] != "COMPLETE":
        raise RuntimeError(result.get("msg", result["status"]))
    table = pd.DataFrame(result["data"], columns=[s["name"] for s in result["fields"]])
    table.to_csv(name, index=False)
    return name, len(table)


if __name__ == "__main__":
    sets = [
        ("mast_cone.csv", "JWST", "NIRCAM/IMAGE"),
        ("miri_fields.csv", "JWST", "MIRI/IMAGE"),
        ("niriss_fields.csv", "JWST", "NIRISS/IMAGE"),
        ("hst_fields.csv", "HST", None),
    ]
    with ThreadPoolExecutor(4) as pool:
        for filename, n in pool.map(query, sets):
            print(filename, n)
