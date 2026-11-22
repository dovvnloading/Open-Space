"""Re-vet only saved untyped NED cones offline, without changing the survey ledgers."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
from unittest.mock import patch

import vet

ONLY_ERROR = "No reported object types in a nonempty NED cone"
ORIGINAL_DECODE = vet.decode_table


def decode(raw, required, expected_count, label=None):
    try:
        return ORIGINAL_DECODE(raw, required, expected_count, label)
    except ValueError as exc:
        # This error can only follow the original successful TAP, schema,
        # canonical-type-field, and row-count checks, with all row types blank.
        if label != "ned" or str(exc) != ONLY_ERROR or "pretype" not in required:
            raise
    columns = tuple(name for name in required if name != "pretype")
    rows = ORIGINAL_DECODE(raw, columns, expected_count)
    for row in rows:
        row["pretype"] = None
    return rows


def replay(root, output):
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != vet.CAMPAIGN_REF
            or os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("Run with the Chandra research workflow")
    root = root.resolve()
    records, errors = [], []
    ledger_paths = sorted((root / "observations").glob("*.json.gz"))
    if not ledger_paths:
        raise ValueError("No deposited observation ledgers")
    with patch.object(vet, "decode_table", decode):
        for path in ledger_paths:
            ledger_bytes = path.read_bytes()
            obs = json.loads(gzip.decompress(ledger_bytes))
            for product in obs.get("products", []):
                search = product.get("search") or {}
                groups = [search.get("triggers", [])] + [c["triggers"] for c in search.get("controls", [])]
                for events in groups:
                    for event in events:
                        original = event.get("vet", {})
                        if original.get("errors", {}).get("ned") != ONLY_ERROR:
                            continue
                        item = {"obsid": obs["obsid"], "trigger_key": event["trigger_key"],
                                "input_sha256": product["input"]["sha256"], "original_status": original["status"],
                                "original_errors": original["errors"], "instrument_pass": event["instrument_pass"],
                                "ledger_sha256": hashlib.sha256(ledger_bytes).hexdigest(),
                                "ledger": str(path.relative_to(root))}
                        records.append(item)
                        try:
                            cache = (root / event.get("catalog_cache", obs["catalog_cache"])).resolve()
                            if not cache.is_relative_to(root):
                                raise ValueError("Saved catalog cache points outside the evidence directory")
                            result = vet.vet_trigger(event, cache, replay=True)
                            for route, receipt in result["queries"].items():
                                if receipt.get("ok"):
                                    old = original["queries"][route]
                                    if (not old.get("ok") or not receipt.get("sha256")
                                            or receipt["sha256"] != old.get("sha256")
                                            or receipt["url"] != old.get("url")):
                                        raise ValueError(f"{route}: catalog input differs from original trigger receipt")
                            item["catalog_cache"] = str(cache.relative_to(root))
                            item["result"] = result
                            if result["errors"]:
                                errors.append([obs["obsid"], event["trigger_key"], result["errors"]])
                        except (KeyError, OSError, ValueError, TypeError) as exc:
                            item["error"] = str(exc)
                            errors.append([obs["obsid"], event["trigger_key"], str(exc)])
    report = {"ledger_count": len(ledger_paths), "affected_triggers": len(records),
              "vet_sha256": hashlib.sha256(Path(vet.__file__).read_bytes()).hexdigest(),
              "replay_implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "interpretation": "Secondary classification from identical catalog response bytes. "
              "Untyped rows participate in the original positional matching as unclassified sources. "
              "Detector decisions, survey FAP and primary ledgers retain their deposited values.",
              "triggers": records, "errors": errors}
    output.parent.mkdir(parents=True, exist_ok=True)
    destination = output.resolve()
    if destination.is_relative_to(root / "observations") or destination.is_relative_to(root / "catalogs"):
        raise ValueError("The replay report needs a separate output path")
    tmp = destination.with_suffix(".tmp")
    tmp.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n")
    tmp.replace(destination)
    print(f"Replayed {len(records)} real or control triggers. Incomplete results {len(errors)}.")
    return bool(errors)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raise SystemExit(replay(args.evidence, args.output))
