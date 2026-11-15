"""Rebuild UTC daily ACE/ULEIS medians and their retrieval receipt."""
import csv
import hashlib
import io
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from statistics import median
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
SERVICE = "https://cdaweb.gsfc.nasa.gov/hapi/"
DATASET = "AC_H2_ULE"
CHANNELS = ("O_S1", "O_L1", "O_L2", "Fe_L1")
QUERY = {
    "id": DATASET, "parameters": "Time," + ",".join(CHANNELS),
    "time.min": "2011-09-18T00:00:00Z", "time.max": "2011-10-04T00:00:00Z",
    "format": "csv",
}


def fetch(url):
    with urlopen(url, timeout=90) as response:
        return response.read()


def analyze(metadata_bytes, data_bytes):
    metadata = json.loads(metadata_bytes)
    instruments = {v["name"]: v for v in metadata["parameters"] if v["name"] in CHANNELS}
    daily = defaultdict(lambda: {key: [] for key in CHANNELS})
    times = set()
    for time, *values in csv.reader(io.StringIO(data_bytes.decode())):
        instant = datetime.fromisoformat(time.replace("Z", "+00:00"))
        if instant in times:
            raise ValueError(f"Duplicate sample at {time}")
        times.add(instant)
        for key, value in zip(CHANNELS, values, strict=True):
            if value and float(value) != float(instruments[key]["fill"]):
                daily[instant.date()][key].append(float(value))
    rows = []
    for date, channel in sorted(daily.items()):
        row = {"date_utc": date.isoformat()}
        for key, values in channel.items():
            row[key + "_median"] = f"{median(values):.8g}" if values else ""
            row[key + "_valid_hours"] = len(values)
        rows.append(row)
    with (ROOT / "daily_medians.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "dataset": DATASET, "data_url": SERVICE + "data?" + urlencode(QUERY),
        "raw_csv_sha256": hashlib.sha256(data_bytes).hexdigest(),
        "raw_rows": len(times), "metadata_url": SERVICE + "info?id=" + DATASET,
        "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
        "channels": {k: {a: instruments[k][a] for a in ["description", "units", "fill"]}
                     for k in CHANNELS},
    }
    (ROOT / "source.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    analyze(fetch(SERVICE + "info?id=" + DATASET), fetch(SERVICE + "data?" + urlencode(QUERY)))
