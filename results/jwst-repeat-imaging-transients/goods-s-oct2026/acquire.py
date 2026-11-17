from s3 import keys, base
from concurrent.futures import ThreadPoolExecutor
import csv, urllib.request as u, pathlib, hashlib
from argparse import ArgumentParser

p = ArgumentParser()
p.add_argument("pid", type=int)
a = p.parse_args()
pid = a.pid
froot = pathlib.Path(str(pid))
froot.mkdir(exist_ok=True)
files = list(keys(f"jwst/public/jw{pid:05d}/L3/t/"))
cat = [(k, n) for k, n in files if k.endswith("_cat.ecsv")]
print(
    "cats:",
    len(cat),
    "totalMB:",
    round(sum(n for k, n in cat) / 1e6, 1),
    "example:",
    cat[:1],
    flush=True,
)


def download(kn):
    k, n = kn
    p = froot / k.split("/")[-1]
    h = hashlib.sha256()
    if not p.exists():
        try:
            data = u.urlopen(base + k, timeout=90).read()
            assert len(data) == n
            p.write_bytes(data)
        except Exception as e:
            print(k.split("/")[-1], str(e), flush=True)
            return
    h.update(p.read_bytes())
    return str(p), base + k, n, h.hexdigest()


with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(download, cat))
with open(f"{pid}_manifest.csv", "w") as f:
    w = csv.writer(f)
    w.writerow(["local_path", "url", "bytes", "sha256"])
    w.writerows(x for x in res if x)
with open(f"{pid}_index.csv", "w") as f:
    w = csv.writer(f)
    w.writerow(["key", "bytes"])
    w.writerows(files)
