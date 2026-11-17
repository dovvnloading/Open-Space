"""Fetch the exact Gaia identifiers in El-Badry et al. (2026), Table 2."""
import csv
import hashlib
import io
import json
import urllib.parse
import urllib.request
import argparse
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--validation", action="store_true", help="Query the full published 92-target sample.")
args = parser.parse_args()
catalogue = "published_targets.csv" if args.validation else "controls.csv"
ids = [r["source_id"] for r in csv.DictReader((HERE / catalogue).open())]
assert len(ids) == len(set(ids)) == (92 if args.validation else 20)
numbers = ", ".join(ids)
queries = {
    "gaia": ("https://gea.esac.esa.int/tap-server/tap/sync", f"""
SELECT source_id, ra, dec, l, b, ecl_lat, parallax, parallax_error,
       pmra, pmdec, pmra_error, pmdec_error, parallax_pmra_corr,
       parallax_pmdec_corr, pmra_pmdec_corr, phot_g_mean_mag, bp_rp,
       nu_eff_used_in_astrometry, pseudocolour, astrometric_params_solved,
       ruwe, visibility_periods_used, ipd_frac_multi_peak,
       ipd_gof_harmonic_amplitude, duplicated_source, has_xp_continuous
FROM gaiadr3.gaia_source
WHERE source_id IN ({numbers})
ORDER BY source_id"""),
    "bj": ("https://dc.g-vo.org/tap/sync", f"""
SELECT source_id, r_med_geo, r_lo_geo, r_hi_geo, r_med_photogeo,
       r_lo_photogeo, r_hi_photogeo
FROM gedr3dist.main
WHERE source_id IN ({numbers})
ORDER BY source_id""")
}
if args.validation:
    # The 92-object comparison needs only Gaia astrometry.
    both = sorted(set(ids) | {r["source_id"] for r in csv.DictReader((HERE / "controls.csv").open())})
    queries = {"validation_gaia": queries["gaia"], "fidelity": ("https://dc.g-vo.org/tap/sync", f"""
SELECT source_id, fidelity_v2, norm_dg, theta_arcsec_worst_source
FROM gedr3spur.main
WHERE source_id IN ({", ".join(both)})
ORDER BY source_id""")}
for stem, (endpoint, query) in queries.items():
    # 200 is safely above the number of distinct IDs.
    query = query.strip() + "\n"
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=200, QUERY=query)
    url = endpoint + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as response:
        data = response.read()
        status, content_type = response.status, response.headers.get("Content-Type")
    table = list(csv.DictReader(io.StringIO(data.decode())))
    returned = [r["source_id"] for r in table]
    requested = set(both if stem == "fidelity" else ids)
    assert len(returned) == len(set(returned)) and set(returned) <= requested
    if stem in ("gaia", "validation_gaia"):
        assert set(returned) == set(ids)
    (HERE / f"{stem}.adql").write_text(query)
    (HERE / f"{stem}.csv").write_bytes(data)
    receipt = dict(endpoint=endpoint, method="GET", url=url, utc=datetime.now(timezone.utc).isoformat(),
                   http_status=status, content_type=content_type, requested_ids=len(requested),
                   rows=len(table), missing_ids=sorted(requested - set(returned)), bytes=len(data),
                   sha256=hashlib.sha256(data).hexdigest(), query_sha256=hashlib.sha256(query.encode()).hexdigest())
    (HERE / f"{stem}_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(stem, len(returned), receipt["sha256"], flush=True)
