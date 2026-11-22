"""Repeat one deposited observation from its exact EVT2 bytes and saved catalog responses."""
import argparse
import math
from pathlib import Path
from tempfile import TemporaryDirectory

from archive import dump, replay
from campaign import EVIDENCE, N_CONTROL, N_INJECT, identity, read, read_obs, response
from search import search_file
from vet import vet_trigger


def reproduce(obsid, temp):
    if identity() != read(EVIDENCE / "inventory/campaign.json")["implementation"]:
        raise ValueError("Use the implementation recorded in the campaign receipt")
    ledger = read_obs(EVIDENCE / "observations" / f"{obsid}.json.gz")
    products = ledger["acquired"]
    dump(temp / "products.json", products)
    replay(temp / "products.json", temp)
    for record in ledger["products"]:
        if not record.get("complete"):
            raise ValueError("This product has no completed reference scan")
        expected = record["search"]
        product = temp / record["input"]["filename"]
        actual = search_file(product, seed=expected["seed"], n_controls=N_CONTROL, n_injections=N_INJECT)
        original_runs = [expected, *expected["controls"]]
        actual_runs = [actual, *actual["controls"]]
        for i, (old, new) in enumerate(zip(original_runs, actual_runs, strict=True)):
            for t, fresh in zip(old["triggers"], new["triggers"], strict=True):
                fresh["obsid"] = obsid
                archived_vet = t.pop("vet")
                cache = t.pop("catalog_cache", ledger["catalog_cache"])
                fresh_vet = vet_trigger(fresh, EVIDENCE / cache, replay=True)
                for key in ("status", "host", "matches"):
                    if fresh_vet[key] != archived_vet.get(key):
                        raise ValueError(f"Counterpart replay differs for {product.name}, {key}")
                t.pop("diagnostics", None)
                flux = t.pop("flux", {})
                key = t.pop("trigger_key")
                if i == 0 and "peak_flux_erg_cm2_s" in flux:
                    values = response(fresh, product, temp / "photometry" / product.name / key)
                    for field in ("peak_flux_erg_cm2_s", "flux_counting_sigma_erg_cm2_s"):
                        if not math.isclose(values[field], flux[field], rel_tol=1e-9):
                            raise ValueError(f"Peak calibration differs for {product.name}, {key}")
        if actual != expected:
            raise ValueError(f"Blind scan, time scrambles or injection replay differs for {product.name}")
    print(f"Reproduced {len(products)} archived products of ObsID {obsid}.")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("obsid", type=int)
    with TemporaryDirectory() as directory:
        reproduce(cli.parse_args().obsid, Path(directory))
