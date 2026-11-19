#!/usr/bin/env python3
"""Independent, read-only replay of the deposited GWTC summary-pair census.

Does NOT import run.py. Independently checks all proxy scores, 1,500
catalogue-stratified shuffles, exact add-one FAP and frozen PE receipts.
NOT a posterior-overlap significance or discovery validation.
"""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

SEED = 541821


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("output_dir", type=Path, nargs="?", default=Path(__file__).resolve().parent / "outputs")
    ap.add_argument("--receipt", type=Path)
    args = ap.parse_args()
    o = args.output_dir
    read = lambda n: json.loads((o/n).read_text())
    pop, bg, manifest, delay = (read(n) for n in (
        "eligible_events.json", "proxy_background.json",
        "release_manifest.json", "time_delay_population.json"))
    with gzip.open(o/"all_pairs_proxy.tsv.gz", "rt") as f:
        recorded = list(csv.DictReader(f, delimiter="\t"))
    names = [x["name"] for x in pop]
    n = len(pop)
    assert len(names) == len(set(names))
    assert len(recorded) == n*(n-1)//2
    bins = np.array(delay["log10_delta_t_bin_edges"])
    pl, pu = (np.array(delay[k]) for k in ("lensed_bin_prob", "unlensed_bin_prob"))
    assert np.isclose(pl.sum(), 1) and np.isclose(pu.sum(), 1)
    assert np.all(pl > 0) and np.all(pu > 0)
    gps = np.array([x["gps"] for x in pop])
    mass, smass, spin, sspin = (np.array([x[k] for x in pop]) for k in (
        "lnmc", "sigma_lnmc", "chieff", "sigma_chieff"))
    intrinsic, D = np.full((n,n), -1e6), np.full((n,n), -1e6)
    calculated = {}
    for i,j in itertools.combinations(range(n),2):
        delta = abs(gps[i]-gps[j])/86400
        idx = int(np.clip(np.digitize(math.log10(max(.001,delta)), bins)-1, 0, len(pl)-1))
        dr = float(math.log(pl[idx]/pu[idx]))
        v, w = smass[i]**2+smass[j]**2, sspin[i]**2+sspin[j]**2
        lm = -.5*((mass[i]-mass[j])**2/v + math.log(2*math.pi*v))
        ls = -.5*((spin[i]-spin[j])**2/w + math.log(2*math.pi*w))
        intrinsic[i,j] = intrinsic[j,i] = lm+ls
        D[i,j] = D[j,i] = dr
        calculated[tuple(sorted((names[i], names[j])))] = (float(delta), float(lm+ls+dr))
    max_error = 0.
    for row in recorded:
        key = tuple(sorted((row["event_1"],row["event_2"])))
        assert key in calculated
        delta, score = calculated.pop(key)
        assert abs(delta-float(row["delay_days"])) < 1e-8
        max_error = max(max_error, abs(score-float(row["proxy_logscore"])))
    assert not calculated and max_error < 1e-10
    upper = np.triu_indices(n,1)
    original_max = float(np.max((intrinsic+D)[upper]))
    assert abs(original_max-float(recorded[0]["proxy_logscore"])) < 1e-10
    assert abs(original_max-bg["max_observed_proxy_logscore"]) < 1e-10
    groups = {}
    for k,e in enumerate(pop): groups.setdefault(e["catalog"], []).append(k)
    rng = np.random.default_rng(SEED+2)
    num = bg["simulated_unlensed_trials"]
    null_max = []
    for _ in range(num):
        perm = np.arange(n)
        for idx in groups.values(): perm[idx] = rng.permutation(idx)
        null_max.append(float(np.max((intrinsic[np.ix_(perm,perm)]+D)[upper])))
    count = sum(x >= original_max for x in null_max)
    p = (count+1)/(num+1)
    assert abs(p-bg["global_proxy_fap"]) < 1e-12
    assert np.allclose(np.quantile(null_max, [.5,.9,.95,.99,1]),
                       bg["null_max_logscore_quantiles"], rtol=1e-12, atol=1e-12)
    receipts = read("posterior_receipts.json") if (o/"posterior_receipts.json").exists() else []
    verified = []
    for r in receipts:
        index = manifest["files"][r["release"]][r["event"]]
        assert index["checksum"].rsplit(":",1)[-1] == r["zenodo_md5"]
        assert index["size"] == r["original_hdf5_size"]
        assert len(r["original_hdf5_sha256"]) == 64
        pth = o/"extracted"/(r["event"]+".npz")
        assert hashlib.sha256(pth.read_bytes()).hexdigest() == r["extracted_npz_sha256"]
        verified.append(r["event"])
    from collections import Counter
    result = {
      "verdict": "PASS: independent all-pair proxy and shuffled FAP exactly reproduced; NOT a global posterior-overlap Bayes-factor FAP",
      "n_eligible": n, "eligible_by_release": dict(Counter(x["release"] for x in pop)),
      "n_unique_pairs": len(recorded), "max_score_absolute_error": max_error,
      "n_null_shuffles": num, "null_exceedances": count, "add_one_proxy_global_fap": p,
      "strongest_pair": [recorded[0]["event_1"],recorded[0]["event_2"]],
      "source_posterior_receipts_verified": len(verified),
      "caveat": "Original HDF5 SHA256 values reported by runner after matching Zenodo MD5; referee verified extracted NPZ SHA256, not redownloaded original HDF5. No sky/spin BF calibration, waveform phase, strain or EM check."
    }
    if args.receipt: args.receipt.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
