from pathlib import Path
import hashlib, urllib.request, concurrent.futures

ROOT = Path(__file__).resolve().parent
SOURCE = "https://www.astrouw.edu.pl/ogle"


def urls():
    for cloud in ["lmc", "smc"]:
        for version, stem in [(3, "ogle3/OIII-CVS"), (4, "ogle4/OCVS")]:
            for item in ["phot.tar.gz", "ident.dat", "remarks.txt", "README"] + (
                ["cep1O2O.dat", "cepF1O.dat"] if version == 4 else []
            ):
                yield f"{cloud}{version}_{item}", f"{SOURCE}/{stem}/{cloud}/cep/{item}"


def get(item):
    name, url = item
    p = ROOT / "inputs" / name
    if not p.exists():
        for i in range(3):
            try:
                with urllib.request.urlopen(url, timeout=90) as response:
                    raw = response.read()
                p.write_bytes(raw)
                break
            except OSError:
                if i == 2:
                    raise
    return (name, hashlib.sha256(p.read_bytes()).hexdigest(), url, p.stat().st_size)


if __name__ == "__main__":
    (ROOT / "inputs").mkdir(exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(6) as pool:
        rows = list(pool.map(get, urls()))
    print("\n".join(" ".join(map(str, r)) for r in rows))
