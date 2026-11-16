"""Retrieve and audit paired RiverSP v2 and SoS v3 time series from PO.DAAC.

Usage: python fetch_swot.py OUTPUT_DIRECTORY REACH_ID [REACH_ID ...]
The output manifest names and hashes the HTTP payloads and the CSV strings
decoded from them. A new directory preserves earlier retrievals.
"""

import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ENDPOINT = "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
FIELDS = (
    "reach_id,river_name,p_lon,p_lat,p_width,p_length,p_dist_out,p_n_nodes,"
    "p_n_ch_max,p_n_ch_mod,p_dam_id,p_low_slp,cycle_id,pass_id,time_str,"
    "range_start_time,range_end_time,granuleUR,crid,sword_version,"
    "collection_shortname,collection_version,reach_q,reach_q_b,xovr_cal_q,"
    "partial_f,obs_frac_n,n_good_nod,dark_frac,ice_clim_f,ice_dyn_f,"
    "layovr_val,node_dist,loc_offset,xtrk_dist,wse,wse_u,wse_r_u,width,"
    "width_u,slope,slope_u,slope_r_u,slope2,slope2_u,slope2_r_u,"
    "sos_consensus_q,sos_hivdi_q,sos_metroman_q,sos_momma_q,sos_sad_q,"
    "sos_sic4dvar_q,sos_lakeflow_q"
)


def acquire(url):
    req = Request(url, headers={"User-Agent": "Open-Space river hydrology research"})
    with urlopen(req, timeout=120) as response:
        payload = response.read()
        meta = {
            "request_url": url,
            "response_url": response.url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "http_date": response.headers.get("Date"),
            "content_type": response.headers.get("Content-Type"),
            "http_status": response.status,
        }
    return payload, meta


def save_bytes(directory, name, data):
    path = directory / name
    path.write_bytes(data)
    return {"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("reaches", nargs="+", type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    params = {"collection_concept_id": "C3905028734-POCLOUD", "page_size": "200"}
    payload, meta = acquire(CMR + "?" + urlencode(params))
    granules = json.loads(payload)["feed"]["entry"]
    africa, = [item for item in granules if item["title"].startswith("af_sword_v16_SOS_results")]
    meta["body"] = save_bytes(args.output, "cmr-l4.json", payload)
    manifest = {
        "geometry_product": "SWOT_L2_HR_RiverSP_2.0",
        "discharge_product": "SWOT_L4_HR_DAWG_SOS_DISCHARGE_V3",
        "discharge_doi": "10.5067/SWOT-SOS-RD3",
        "africa_results_granule": {
            k: africa[k] for k in ("id", "title", "time_start", "time_end", "updated")
        },
        "cmr": meta,
        "reach_requests": [],
    }
    for reach in args.reaches:
        params = {
            "collection_name": "SWOT_L2_HR_RiverSP_2.0",
            "feature": "Reach",
            "feature_id": str(reach),
            "start_time": "2023-03-01T00:00:00Z",
            "end_time": "2025-05-04T00:00:00Z",
            "output": "csv",
            "fields": FIELDS,
        }
        payload, meta = acquire(ENDPOINT + "?" + urlencode(params))
        response = json.loads(payload)
        if response.get("status") != "200 OK" or not response.get("hits"):
            raise ValueError(f"Reach {reach}: {response}")
        csv_bytes = response["results"]["csv"].encode("utf-8")
        records = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8"))))
        if len(records) != response["hits"] or any(row["reach_id"] != str(reach) for row in records):
            raise ValueError(f"Row count or reach identifier mismatch for {reach}")
        missing_fields = set(FIELDS.split(",")) - set(records[0])
        if missing_fields:
            raise ValueError(f"Missing fields for {reach}: {missing_fields}")
        meta.update({
            "reach": str(reach),
            "hits": response["hits"],
            "body": save_bytes(args.output, f"{reach}.json", payload),
            "decoded_csv": save_bytes(args.output, f"{reach}.csv", csv_bytes),
        })
        manifest["reach_requests"].append(meta)
        print(f"{reach}: {len(records)} observations", flush=True)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
