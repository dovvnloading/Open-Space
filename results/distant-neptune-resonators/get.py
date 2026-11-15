from pathlib import Path
import json, urllib.request, urllib.parse, time

root = Path(__file__).parent
store = root / "work"


def save(url, name):
    dest = store / name
    if not dest.exists():
        for n in range(4):
            try:
                with urllib.request.urlopen(url, timeout=90) as response:
                    data = response.read()
                if name.endswith(".json"):
                    json.loads(data)
                dest.write_bytes(data)
                break
            except Exception:
                if n == 3:
                    raise
                time.sleep(2)
    print(name, dest.stat().st_size, flush=True)
    return dest


if __name__ == "__main__":
    q = {
        "sb-kind": "a",
        "fields": "spkid,full_name,pdes,name,a,e,i,om,w,ma,q,epoch,data_arc,n_obs_used",
        "sb-cdata": json.dumps(
            {"AND": ["a|RG|113|400", "q|RG|29|60", "data_arc|GT|720"]}
        ),
        "sort": "a",
        "full-prec": "true",
    }
    save(
        "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urllib.parse.urlencode(q),
        "catalog.json",
    )
    save(
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
        "de440s.bsp",
    )
