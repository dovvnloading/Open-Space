"""Fetch the fast branch. CSV byte hashes and the exact ADQL are retained beside each response."""
import csv
import hashlib
import json
import pathlib
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parent
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
GAVO = "https://dc.g-vo.org/tap/sync"
SIMBAD = "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync"


def fetch(name, endpoint, query):
    query = query.strip()
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)
    url = endpoint + "?" + urllib.parse.urlencode(params)
    started = datetime.now(timezone.utc).isoformat()
    with urllib.request.urlopen(url, timeout=180) as response:
        payload = response.read()
        content_type = response.getheader("Content-Type", "").lower()
        if not payload or payload.lstrip().startswith(b"<") or "csv" not in content_type:
            raise ValueError(f"{name}: invalid CSV response ({content_type}, {payload[:160]!r})")
        rows = list(csv.DictReader(payload.decode().splitlines()))
        receipt = dict(endpoint=endpoint, method="GET", query=query, params=params, status=response.status,
                       server_date=response.getheader("Date"), content_type=content_type,
                       started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
                       rows=len(rows), nbytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    (ROOT / (name + ".adql")).write_text(query + "\n")
    (ROOT / (name + ".csv")).write_bytes(payload)
    (ROOT / (name + "_receipt.json")).write_text(json.dumps(receipt, indent=2) + "\n")
    return rows


def followup(rows):
    ids = ",".join(row["source_id"] for row in rows)
    names = ",".join("'Gaia DR3 " + row["source_id"] + "'" for row in rows)
    cols = """source_id, ra, dec, parallax, parallax_error, parallax_over_error, pm,
pmra, pmra_error, pmdec, pmdec_error, parallax_pmra_corr, parallax_pmdec_corr,
pmra_pmdec_corr, phot_g_mean_mag, bp_rp, ruwe, duplicated_source,
visibility_periods_used, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude,
astrometric_params_solved, nu_eff_used_in_astrometry, pseudocolour, ecl_lat,
astrometric_excess_noise_sig, phot_bp_rp_excess_factor, has_xp_continuous, radial_velocity"""
    followups = [
        ("gaia", GAIA, f"SELECT {cols} FROM gaiadr3.gaia_source WHERE source_id IN ({ids})"),
        ("distance", GAVO, f"SELECT source_id,r_med_geo,r_lo_geo,r_hi_geo,r_med_photogeo,"
         f"r_lo_photogeo,r_hi_photogeo FROM gedr3dist.main WHERE source_id IN ({ids})"),
        ("fidelity", GAVO, f"SELECT source_id,fidelity_v2,theta_arcsec_worst_source,norm_dg,"
         f"dist_nearest_neighbor_at_least_0_brighter FROM gedr3spur.main WHERE source_id IN ({ids})"),
        ("simbad", SIMBAD, f"SELECT i.id, b.main_id, b.otype, b.sp_type, b.rvz_radvel, b.nbref "
         f"FROM ident AS i JOIN basic AS b ON i.oidref=b.oid WHERE i.id IN ({names})"),
    ]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = [pool.submit(fetch, *query) for query in followups]
        for i, future in enumerate(results):
            data = future.result()
            if i < 3 and {r["source_id"] for r in data} != set(ids.split(",")):
                raise RuntimeError(f"{followups[i][0]} returned a different ID set")


def main():
    where = (ROOT / "count.adql").read_text().strip().split("\n", 1)[1]
    query = "SELECT source_id,ra,dec,parallax,parallax_error,pm,ruwe,bp_rp,phot_g_mean_mag\n" + where
    with ThreadPoolExecutor(max_workers=2) as pool:
        counts = pool.submit(fetch, "count", GAIA, "SELECT COUNT(*) AS n\n" + where)
        parent = pool.submit(fetch, "parent", GAIA, query)
        n = int(counts.result()[0]["n"])
        rows = parent.result()
    if len(rows) != n or len({r["source_id"] for r in rows}) != n:
        raise RuntimeError("Parent enumeration and separate COUNT disagree")
    followup(rows)
    # Replay the frozen coordinates and IDs used for ancillary checks.
    requests = [("gaia_3arcsec", GAIA), ("dr2", GAIA), ("simbad_5arcsec", SIMBAD),
                ("ap", GAIA), ("negative_count", GAIA), ("templates", GAIA), ("eligible_count", GAIA)]
    with ThreadPoolExecutor(max_workers=3) as pool:
        pending = [pool.submit(fetch, name, endpoint, (ROOT / (name + ".adql")).read_text())
                   for name, endpoint in requests]
        for future in pending:
            future.result()


if __name__ == "__main__":
    main()
