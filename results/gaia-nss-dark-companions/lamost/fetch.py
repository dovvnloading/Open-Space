"""Fetch DR11 v2.0 MRS exact DR3-ID matches through the public catalogue form."""
import argparse
import csv
import gzip
import hashlib
import http.cookiejar
import io
import json
import re
import time
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

URL = "https://www.lamost.org/dr11/v2.0/table/med_catalogue"
COLUMNS = """mobsid obsid designation obsdate mjd planid spid fiberid lmjm band ra_obs dec_obs snr gaia_source_id fibertype offsets rv_b1 rv_b1_err rv_b_flag rv_r1 rv_r1_err rv_r_flag rv_br1 rv_br1_err rv_br_flag coadd fibermask bad_b bad_r moon_flg subproject""".split()
CONTROLS = {"4373465352415301632", "5870569352746779008", "4318465066420528000", "6328149636482597888"}


def query(ids, table="med_catalogue", columns=COLUMNS):
    form_url = URL.rsplit("/",1)[0] + "/" + table
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    with opener.open(form_url, timeout=70) as page:
        token = re.search(rb'name="csrf_token" value="([^"]+)"', page.read()).group(1).decode()
    form = [
        ("csrf_token", token), ("pos_group", "ra,dec"),
        (table + ".gaia_source_id.operation", "in"),
        (table + ".gaia_source_id.textarea", "\n".join(ids)),
    ] + [("showcol", table + "." + col) for col in columns]
    req = urllib.request.Request(form_url + "/q", urllib.parse.urlencode(form).encode(), headers={"Referer": form_url})
    with opener.open(req, timeout=100) as page:
        path = re.search(rb'var save_as_base_url = "([^"]+)"', page.read()).group(1).decode()
    url = urllib.parse.urljoin(form_url, path) + "output.fmt=csv&"
    with opener.open(url, timeout=180) as response:
        raw = response.read()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if set(columns) - set(reader.fieldnames or []):
        raise ValueError("Expected MRS CSV header missing")
    rows = list(reader)
    if any(row["gaia_source_id"] not in ids for row in rows):
        raise ValueError("Service returned unrequested ID")
    return url, raw, rows


def bundle(directory, receipts):
    """Collect exact response files into one ZIP without changing their bytes."""
    destination = directory / "responses.zip"
    temporary = destination.with_suffix(".part")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_STORED) as out:
        for receipt in receipts:
            name = receipt["file"]
            out.writestr(zipfile.ZipInfo(name, date_time=(1980,1,1,0,0,0)),
                         (directory / name).read_bytes())
    temporary.replace(destination)
    for receipt in receipts:
        (directory / receipt["file"]).unlink()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mass_screen", type=Path)
    p.add_argument("directory", type=Path)
    p.add_argument("--bundle", action="store_true", help="Store exact responses in one ZIP")
    a = p.parse_args()
    a.directory.mkdir(parents=True, exist_ok=True)
    with a.mass_screen.open() as f:
        ids = sorted({row["source_id"] for row in csv.DictReader(f)} | CONTROLS, key=int)
    groups = [ids[i:i+100] for i in range(0, len(ids), 100)]

    def process(item):
        n, group = item
        for attempt in range(3):
            try:
                url, data, rows = query(group)
                break
            except (ValueError, OSError, AttributeError):
                if attempt == 2:
                    raise
                time.sleep(3 * (attempt + 1))
        name = f"response-{n:02d}.csv.gz"
        (a.directory / name).write_bytes(gzip.compress(data, mtime=0))
        return dict(file=name, url=url, sha256_csv=hashlib.sha256(data).hexdigest(),
                    requested_ids=group, rows=len(rows), matched_ids=len({r['gaia_source_id'] for r in rows}))

    with ThreadPoolExecutor(max_workers=3) as pool:
        receipts = list(pool.map(process, enumerate(groups)))
    (a.directory / "receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")
    if a.bundle:
        bundle(a.directory, receipts)


if __name__ == "__main__":
    main()
