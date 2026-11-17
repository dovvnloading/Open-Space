# A/2024 U2 is bound on entry and two marginal CODE comets admit near-parabolic solutions

We independently fitted 337 published ground positions for three objects with small nominal excess speeds. The asteroid-like A/2024 U2 was observed over 721 days. The two shorter comet arcs were selected because they have negative published original reciprocal semimajor axes. Both the observed coordinates and their dates come from the November 17, 2026 live MPC API queries.

| Object | Observations, used/received | Gravity-fit RMS per coordinate (arcsec) | Gravity-only inbound 1/a at 250 au (au⁻¹, formal 1σ) | Measured sensitivity |
|---|---:|---:|---:|---|
| A/2024 U2 | 159/159 | 0.445 | +0.00002521 ± 0.00000395 | Bound at 250 au despite a nominally hyperbolic heliocentric catalogue orbit. The current arc ends 2026 Oct 17. |
| C/2011 M1 (LINEAR) | 117/122 | 1.208 | −0.002036 ± 0.000703 | The outward Marsden radial/transverse fit gives −0.001284 ± 0.002711. Forcing a parabolic incoming orbit costs Δχ² = 0.43 under that model. |
| C/2020 W5 (Lemmon) | 56/56 | 0.519 | −0.0001957 ± 0.0001486 | A parabolic fit costs Δχ² = 1.66 under gravity and 0.13 with a free radial inverse-square force. |

Positive 1/a denotes a bound orbit. All quoted original values use an actual inbound crossing of 250 au with negative radial velocity. The profile tests constrain 1/a to zero at the inbound 50-au sphere, matching the method used for C/2025 V1.

Removing the final night of C/2011 M1 **reverses the gravity-only energy classification**. Its 250-au result becomes +0.002315 ± 0.001102 au⁻¹, and the held-out positions have a 6.16″ RMS per coordinate. Removing the first night of C/2020 W5 gives −0.0000035 ± 0.0002359 au⁻¹. The formal Marsden coefficients for W5 extrapolate a water-ice activity law from a narrow interval near 3.36 au and have extreme amplitudes at 1 au. They have no demonstrated physical interpretation.

CODE's gravity solutions give −0.001547 ± 0.000506 au⁻¹ (M1) and −0.000203 ± 0.000059 au⁻¹ (W5). Their nominal values agree with the independent gravity-only estimates to within the quoted uncertainties. The catalogue and this analysis use different astrometric weights and rejection procedures. [Published original-orbit crossmatch](../../catalog/code-original.csv) · [C/2011 M1 in CODE](https://pad2.astro.amu.edu.pl/CODE/orbit.php?int=2011m1a5&orb=original) · [C/2020 W5 in CODE](https://pad2.astro.amu.edu.pl/CODE/orbit.php?int=2020w5a2&orb=original)

## Reproduce

Run `python replay.py` from this directory with the dependencies listed in the [independent fitter documentation](../README.md). The script fetches the selected MPC observations and starting vectors, fits the sky directions and propagates the results. The [frozen results and receipts](evidence/results.json) record the query URLs and payloads, wire and observation SHA-256 hashes, fitted states, covariances, residual summaries, holdouts and energy profiles. The three compressed observation files preserve the exact measured positions and times.

To replay those measurements after the live catalogue changes, populate the local observation-response cache from the frozen records before invoking the script. The ordinary fit uses a fresh Horizons initialization if a starting vector has not already been cached. All six state components are free during fitting.

```sh
python - <<'PY'
import gzip, hashlib, json
from pathlib import Path
cache = Path("cache")
cache.mkdir(exist_ok=True)
reference = json.loads(Path("evidence/results.json").read_text())
for obj, row in reference["objects"].items():
    raw = gzip.decompress(Path("evidence", obj + ".obs.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == row["input"]["sha256"]
    (cache / (obj + ".json")).write_text(json.dumps([{"OBS80": raw.decode("ascii")}]))
PY
python replay.py
```
