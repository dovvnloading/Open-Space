"""Capture unmodified Hydrocron response and reproducible manifest for selected Congo reaches."""

import csv, concurrent.futures, hashlib, io, json, pathlib, time, urllib.parse, urllib.request, datetime

ROOT = pathlib.Path(__file__).resolve().parent
REACHES = [
    "13217000061",
    "13217000091",
    "13217000101",
    "13217000111",
    "13217000331",
    "13217000341",
    "13217000351",
    "13217000391",
    "13217000401",
    "13219000031",
    "13219000041",
]
FIELDS = "reach_id,time_str,wse,wse_u,slope,slope_u,slope2,slope2_u,width,width_u,sos_consensus_q,sos_hivdi_q,sos_metroman_q,sos_momma_q,sos_sad_q,sos_sic4dvar_q,reach_q,reach_q_b,partial_f,n_good_nod,obs_frac_n,loc_offset,xtrk_dist,cycle_id,pass_id,crid,sword_version,granuleUR,collection_shortname,collection_version"


def fetch(rid):
    url = (
        "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?"
        + urllib.parse.urlencode(
            {
                "feature": "Reach",
                "feature_id": rid,
                "output": "csv",
                "start_time": "2023-03-29T00:00:00Z",
                "end_time": "2025-09-01T00:00:00Z",
                "collection_name": "SWOT_L2_HR_RiverSP_2.0",
                "fields": FIELDS,
            }
        )
    )
    target = ROOT / "raw" / f"{rid}.json"
    target.parent.mkdir(exist_ok=True)
    if not target.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=45) as response:
                    data = response.read()
                    assert response.status == 200
                parsed = json.loads(data)
                assert parsed["status"] == "200 OK", parsed
                target.write_bytes(data)
                break
            except Exception as exc:
                if attempt == 2:
                    raise
                print(rid, str(exc)[:160], flush=True)
                time.sleep(3)
    data = target.read_bytes()
    rows = list(csv.DictReader(io.StringIO(json.loads(data)["results"]["csv"])))
    print(
        rid, len(rows), sum(float(r["sos_consensus_q"]) >= 0 for r in rows), flush=True
    )
    return {
        "filename": "raw/" + target.name,
        "query_url": url,
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "retrieved_utc": datetime.datetime.fromtimestamp(
            target.stat().st_mtime, datetime.timezone.utc
        ).isoformat(),
        "rows": len(rows),
    }


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(5) as pool:
        entries = list(pool.map(fetch, REACHES))
    (ROOT / "manifest.json").write_text(json.dumps(entries, indent=2) + "\n")
