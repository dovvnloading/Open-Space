"""Archive fixed-node geometry to test whether changing reach support drives the head contrast."""

import concurrent.futures, csv, datetime, hashlib, io, json, pathlib, urllib.parse, urllib.request

ROOT = pathlib.Path(__file__).resolve().parent
# Central fixed positions and flanking positions within each reach, indexed downstream to upstream.
NODES = {"13217000061": [15, 30, 45], "13217000091": [12, 23, 35]}
fields = "reach_id,node_id,time_str,wse,wse_u,lat,lon,node_q,node_q_b,xtrk_dist,cycle_id,pass_id,crid,sword_version,granuleUR"


def fetch(node):
    url = (
        "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?"
        + urllib.parse.urlencode(
            dict(
                feature="Node",
                feature_id=node,
                output="csv",
                start_time="2023-03-29T00:00:00Z",
                end_time="2025-09-01T00:00:00Z",
                collection_name="SWOT_L2_HR_RiverSP_2.0",
                fields=fields,
            )
        )
    )
    target = ROOT / "raw" / f"node-{node}.json"
    if not target.exists():
        with urllib.request.urlopen(url, timeout=50) as response:
            data = response.read()
        assert json.loads(data)["status"] == "200 OK"
        target.write_bytes(data)
    data = target.read_bytes()
    print(node, json.loads(data)["hits"], flush=True)
    return {
        "filename": "raw/" + target.name,
        "query_url": url,
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "retrieved_utc": datetime.datetime.fromtimestamp(
            target.stat().st_mtime, datetime.timezone.utc
        ).isoformat(),
        "rows": json.loads(data)["hits"],
    }


if __name__ == "__main__":
    ids = [rid[:-1] + f"{v:03d}" + rid[-1] for rid, ns in NODES.items() for v in ns]
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        rows = list(pool.map(fetch, ids))
    (ROOT / "node-manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
