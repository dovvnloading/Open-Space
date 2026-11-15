"""Replay the MPC alias/astrometry and ESA orbit checks. Uses the standard library."""
import json
import math
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor


PAIRS = [
    ("2026 SW11", "2026 SA12"),
    ("2026 VO", "2026 VR1"),
    ("2017 RB15", "2026 RG39"),
    ("2007 VE138", "2026 VV15"),
]
CONTROL = ("2018 NC15", "2026 LD1")


def read(url, body=None):
    if body is not None:
        request = urllib.request.Request(
            url, data=json.dumps(body).encode(), method="GET",
            headers={"Content-Type": "application/json"},
        )
    else:
        request = url
    with urllib.request.urlopen(request, timeout=40) as response:
        return response.read()


def query(endpoint, body):
    return json.loads(read("https://data.minorplanetcenter.net/api/" + endpoint, body))


def observations(name):
    result = query("get-obs", {"desigs": [name], "output_format": ["OBS80"]})
    return result[0]["OBS80"].splitlines()


def sky(line):
    h, m, sec = map(float, line[32:44].split())
    d, dm, ds = map(float, line[44:56].split())
    return (
        math.radians(15 * (h + m / 60 + sec / 3600)),
        math.radians(math.copysign(abs(d) + dm / 60 + ds / 3600, d)),
    )


def separation(a, b):
    ra1, dec1 = sky(a)
    ra2, dec2 = sky(b)
    u = math.sin((dec2 - dec1) / 2) ** 2
    u += math.cos(dec1) * math.cos(dec2) * math.sin((ra2 - ra1) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(min(u, 1.0)))) * 3600


def coincident(a, b):
    # The time and station are independently stored in both OBS80 records.
    samples = defaultdict(list)
    for line in a:
        samples[(line[15:32].strip(), line[77:80])].append(line)
    matches = []
    for line in b:
        instant = line[15:32].strip(), line[77:80]
        for earlier in samples[instant]:
            matches.append({"date_utc": instant[0], "station": instant[1],
                            "separation_arcsec": round(separation(earlier, line), 3)})
    return matches


def esa_orbit_name(name):
    url = "https://neo.ssa.esa.int/PSDB-portlet/download?file=" + name.replace(" ", "") + ".ke0"
    body = read(url).decode()
    if not any(row.startswith(" KEP ") for row in body.splitlines()):
        raise ValueError("No orbital elements returned for " + name)
    return body.split("END_OF_HEADER", 1)[1].strip().splitlines()[0]


def main():
    names = [name for pair in PAIRS for name in pair]
    identity = query("query-identifier", {"ids": names})
    all_names = names + list(CONTROL)
    with ThreadPoolExecutor(max_workers=6) as pool:
        obs = dict(zip(all_names, pool.map(observations, all_names)))
        orbits = dict(zip(names, pool.map(esa_orbit_name, names)))
    print(json.dumps({
        "official_records": {
            name: {
                "mpc_primary": identity[name]["unpacked_primary_provisional_designation"],
                "mpc_secondary": identity[name]["unpacked_secondary_provisional_designations"],
                "esa_orbit_file_object": orbits[name],
                "mpc_observations": len(obs[name]),
            } for name in names
        },
        "same_station_same_time": {
            a + " / " + b: coincident(obs[a], obs[b]) for a, b in [*PAIRS, CONTROL]
        },
    }, indent=2))


if __name__ == "__main__":
    main()
