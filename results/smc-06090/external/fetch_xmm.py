"""Retrieve the exact XMM-Newton OM source lists used in the archived comparison."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent

def fetch():
    for record in json.loads((HERE / "xmm.json").read_text())["files"]:
        target = HERE / "inputs" / record["name"]
        if target.exists():
            data = target.read_bytes()
        else:
            with urlopen(record["url"], timeout=90) as response:
                data = response.read()
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError(f"Unexpected SHA-256 for {target.name}")
        if not target.exists():
            target.write_bytes(data)
        print(target.name)

if __name__ == "__main__":
    fetch()
