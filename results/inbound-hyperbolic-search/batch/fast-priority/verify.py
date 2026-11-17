"""Verify archived observations, selected catalogue IDs, source hashes and the roving cross-check."""
import gzip
import json
from pathlib import Path

from replay import KNOWN, cat, hashbytes, jd, vinf

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"


def verify_group(prefix, parts, expected):
    records = {}
    for part in range(parts):
        key = f"{prefix}shard-{part}"
        run = json.loads((EVIDENCE / (key + ".json")).read_text())
        bundle = EVIDENCE / f"{prefix}observations-{part}.json.gz"
        inputs = json.loads(gzip.decompress(bundle.read_bytes()))
        receipt = run["receipt"]
        assert (receipt["shard"], receipt["parts"]) == (part, parts)
        assert receipt["sha256_compressed"] == "680c9ee88ef6f615c30ca4faa347a77d5bed5787c8bbc74d8c7257c62a2aa54f"
        assert receipt["kernel_sha256"] == "dc8fa8136e71df25c1e4183cd777e6ae536eac099dc6225ad77192cf2192221d"
        assert receipt["fitter_sha256"] == "7c126536dafd0c8f42008ba26a6d143624f6626993e0bd523c67603b592699ba"
        assert set(inputs) == {r["spkid"] for r in run["records"]}
        for result in run["records"]:
            sk = result["spkid"]
            assert sk not in records
            raw = inputs[sk].encode("ascii")
            assert hashbytes(raw) == result["raw_sha256"]
            assert len(raw.splitlines()) == result["input_lines"]
            if "staged_sha256" in result:
                header = b"#interstellar\n"
                if prefix == "control-":
                    header += b"#Combine all\n"
                assert hashbytes(header + raw) == result["staged_sha256"]
            if result.get("far_inbound_100au"):
                assert result["status"] == "identity and observation epochs verified"
                assert result["center"] == "Object_-1" and result["inbound"] and result["r_au"] >= 100
                assert abs(result["requested_epoch"] - result["returned_epoch"]) <= 0.0001
                assert result["far_inbound_250au"] == (result["r_au"] >= 250)
            records[sk] = result
    assert set(records) == expected
    return records


def verify_rover(review_file, raw_file, results, n_pairs):
    review = json.loads((EVIDENCE / review_file).read_text())
    result = results[review["spkid"]]
    raw = json.loads(gzip.decompress((EVIDENCE / raw_file).read_bytes()))[review["spkid"]]
    assert review["input_sha256"] == result["raw_sha256"]
    assert result["station_time_mismatches"] == len(review["position_pairs"]) == n_pairs
    sites = list(dict.fromkeys(line[32:65] for line in raw.splitlines() if line[14] == "v"))
    assert len(sites) == 2 and review["first_rover_site"][32:65] == sites[0]
    pairs = {(site, t) for site, t in result["mismatch_examples"]}
    for pair in review["position_pairs"]:
        obs, pos = pair["position"], pair["site_coordinates"]
        assert obs in raw.splitlines() and pos in raw.splitlines()
        assert obs[14] == "V" and pos[14] == "v" and pos[32:65] == sites[1]
        assert obs[77:80] == pos[77:80] == "247" and pair["returned_station"] == "24a"
        assert (pair["returned_station"], pair["returned_utc_jd"]) in pairs
        assert jd(obs.encode("ascii")) == jd(pos.encode("ascii"))
        assert abs(jd(obs.encode("ascii")) - pair["returned_utc_jd"]) < 0.000001


def main():
    catalog, _, catalog_sha = cat()
    assert catalog_sha == "90e391ede544a92c4648436656f38a42cc81f7d2ef7fea3b70dd832c24a26022"
    expected = {str(r["spkid"]) for r in catalog if vinf(r["e"], r["q"]) >= 1}
    assert len(expected) == 67
    primary = verify_group("", 6, expected)
    names = set((HERE / "soho-targets.txt").read_text().splitlines())
    supplemental = {str(r["spkid"]) for r in catalog if r["pdes"] in names}
    assert len(supplemental) == 4 and not (expected & supplemental)
    verify_group("supplement-", 2, supplemental)
    controls = verify_group("control-", 3, set(KNOWN))
    for sk, result in controls.items():
        assert result["raw_sha256"] == primary[sk]["raw_sha256"]

    verify_rover("station-review.json", "observations-3.json.gz", primary, 2)
    verify_rover("control-station-review.json", "control-observations-1.json.gz", controls, 3)
    print("Verified 67 priority objects, four disjoint SOHO records, three ISO alias replays and five rover pairs.")


if __name__ == "__main__":
    main()
