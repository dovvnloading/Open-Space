#!/usr/bin/env python3
"""Independent standard-library checks of the frozen *scope*, not raw data."""
import datetime as dt
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
m=json.loads((HERE/"frozen_sources.json").read_text())
s=m["dja"]["spectral_shards"]
assert m["schema"]=="open-space/jwst-nirspec-frozen-source-v1"
assert m["dja"]["version"]=="v4.4"
assert m["dja"]["expected_entries"]==80367
assert len(s)==9 and sum(item["expected_spectra"] for item in s)==80367
assert len({x["band"] for x in s})==9
assert sum(x["band"]=="prism" for x in s)==1
for item in s+[m["dja"]["catalog"],m["dja"]["columns"]]:
    assert re.fullmatch("[0-9a-f]{32}",item["publisher_md5"])
    assert item["url"].startswith("https://zenodo.org/records/15472354/files/")
    assert item["verified"] is False
    assert item["local_sha256"] is None
rows=m["mast"]["rows"]
assert len(rows)==16
assert len({(r["program"],r["target"]) for r in rows})==len(rows)
assert all(dt.date.fromisoformat(r["release_as_displayed"][:10]) <=
           dt.date.fromisoformat(m["frozen_on"]) for r in rows)
assert all(r["public_product_query_status"]=="not_executed" for r in rows)
assert m["dja"]["download_attempt"]["verified_sha256_count"]==0
assert m["claim_status"].startswith("UNTESTED")
assert not (HERE/"screen_rows.csv").exists()
d={"dja_release":m["dja"]["version"],"expected_spectra":80367,
   "shards":len(s),"mast_dashboard_rows":len(rows),
   "mast_dashboard_reported_count":sum(r["reported_count"] for r in rows),
   "raw_FITS_verified":0,"actual_redshift_fits":0,
   "finding":"PRE-FLIGHT ONLY; NO SURVEY-WIDE NULL OR DISCOVERY"}
print(json.dumps(d,sort_keys=True))
