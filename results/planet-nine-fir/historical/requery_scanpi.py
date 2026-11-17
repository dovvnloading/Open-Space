"""Submit one archived IRSA SCANPI request and obtain a fresh results page."""
import argparse
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("request", type=Path, help="scanpi/request.json or scanpi/galaxy_request.json")
    args = parser.parse_args()
    record = json.loads(args.request.read_text())
    boundary = "ScanpiReproduction"
    fields = []
    for key, value in record["form"].items():
        fields.extend([f"--{boundary}", f'Content-Disposition: form-data; name="{key}"', "", value])
    body = "\r\n".join(fields + [f"--{boundary}--", ""]).encode()
    submit = Request(record["url"], data=body,
                     headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urlopen(submit, timeout=120) as result:
        monitor = result.url
        page = result.read().decode()
    print("Monitor this submitted job at", monitor, flush=True)
    for attempt in range(80):
        match = re.search(r"(https://[^'\"<>]+?/summary\.html)", page)
        if match:
            print("Summary and profile links at", match.group(1))
            return
        time.sleep(15)
        with urlopen(monitor, timeout=120) as result:
            page = result.read().decode()
    raise RuntimeError(f"Job still pending. Resume monitoring {monitor}.")


if __name__ == "__main__":
    main()
