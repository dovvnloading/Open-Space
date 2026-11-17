# High-speed catalogue replay

We independently screened all 67 objects in the frozen JPL sample with nominal osculating heliocentric excess speed at least 1 km/s. The 120-second gravity-only fits produced strong apparent excesses for comet fragments, a nine-day-arc comet and known interstellar objects. None validates a new interstellar object. Four separately chosen SOHO sungrazers reproduce extreme nominal speeds with large astrometric residuals. Their four fits are outside the 67.

## Completed

| Main replay with original automatic alias handling | Objects |
|---|---:|
| Frozen catalogue objects above threshold | 67 |
| At least three associated real rows, sent to Find_Orb | 62 |
| Strict packed-ID and station/time validation | 59 |
| Additional fit verified by matching both internally renamed roving observations | 1 |
| More than one completed output, 2I/Borisov | 1 |
| CPU timeout, C/1898 M1 | 1 |
| No associated rows in the bulk file or the per-designation MPC lookup | 5 |
| Strictly verified barycentric solutions inbound at ≥100 / ≥250 au | 50 / 31 |
| Same boundary counts including the verified rover | 51 / 32 |

The two boundary counts use the estimated radius at the reported epoch. The 50 strictly verified solutions beyond 100 au include 29 with positive nominal energy, 25 with positive nominal energy and a finite positive eccentricity uncertainty, and four with `(e−1)/σe ≥ 5`. All four also reached 250 au. The separately reviewed rover adds one positive ≥5σ fragment at both boundaries. These are local formal errors under one gravitational fit. The universal search supplies no statistically calibrated false-alarm probability. [All 67 records](priority-index.csv) · [Rover epochs, site coordinates and source-code reference](evidence/station-review.json).

| Formal positive outcomes of interest | Excess km/s | (e−1)/σe | Used / returned observations | Used arc, d | Coordinate RMS, arcsec |
|---|---:|---:|---:|---:|---:|
| 1I/'Oumuamua, original automatic alias run | 26.175 | 172.79 | 97 / 105 | 15.75 | 0.441 |
| C/2025 V1 | 3.827 | 6.70 | 87 / 117 | 9.62 | 0.527 |
| C/2025 K1-D | 2.933 | 28.29 | 212 / 227 | 25.74 | 0.413 |
| C/2019 Y4-A | 2.855 | 13.33 | 218 / 348 | 24.34 | 0.665 |
| C/2025 K1-B, rover cross-check | 2.145 | 33.86 | 291 / 358 | 62.67 | 0.627 |
| C/2011 M1 | 1.220 | 3.32 | 76 / 122 | 40.05 | 0.594 |

A post-breakup component has no independent pre-breakup orbit to recover by backward integration of its post-breakup astrometry alone. Hubble imaged the disintegration of both parent comets. The frozen JPL row for C/2025 K1 contains nonzero A1, A2 and A3. The independent work on C/2025 V1 uses its short arc to examine model and observatory dependence. [C/2019 Y4 HST study](https://arxiv.org/abs/2105.02269) · [C/2025 K1 HST study](https://arxiv.org/abs/2511.19707) · [Independent orbit investigations](../../orbits/).

### Known ISO alias checks

The initial automatic run split Borisov's permanent and provisional designations into separate element files. It retained an 8.49-day subset on 3I and a 15.75-day subset on 1I. A second run checked every nonblank permanent and provisional ID in the three committed real-observation bundles. Each bundle contained only that control's published aliases. The native `#Combine all` directive then allowed the fitter to process these accepted aliases together. Each residual still had to match an accepted packed ID, observatory and input epoch. These three repeat fits share the 67-object denominator. [Alias-joined control run](https://github.com/dovvnloading/Open-Space/actions/runs/44037636030).

| Alias-joined control | Excess km/s | Used / returned | Used / input span, days | Coordinate RMS, arcsec | Inbound radius, au |
|---|---:|---:|---:|---:|---:|
| 1I/'Oumuamua | 26.363 | 167 / 220 | 34.51 / 80.04 | 0.326 | 1,676.9 |
| 2I/Borisov | 32.284 | 2,424 / 3,361 | 376.02 / 656.79 | 0.592 | 2,048.4 |
| 3I/ATLAS | 57.968 | 4,100 / 8,151 | 328.37 / 333.44 | 0.599 | 3,670.5 |

The alias-joined 3I output has three literal 247→24a station mismatches in its machine record. The accompanying review matches all three exact epochs and V/v coordinate pairs to the second distinct roving site. Including this check, the native joined-control fits recovered all three known ISO speeds on arcs with substantially different observation counts. After applying the verified rover correction to K1-B and the control-specific replays, 61 distinct targets have an identity-checked fit: 52 inbound beyond 100 au and 33 beyond 250 au. [3I rover rows](evidence/control-station-review.json) · [2I](evidence/control-shard-0.json) · [3I](evidence/control-shard-1.json) · [1I](evidence/control-shard-2.json).

The raw MPC file may include observatory-position continuation rows, duplicate measurements and aliases. Returned observations are the fitter's list before its automatic exclusions. The used count and used time span are both reported. Matching all returned identities and epochs does not establish that the fitter used every input position. The known-control comparison measures this limitation directly. Cases with insufficient astrometry were C/1853 R1, C/1940 R2, C/1940 S1, C/1959 O1 and C/2020 P4-A. C/1898 M1 reached its CPU limit with 27 input lines.

### Four follow-ups from the full catalogue

Each object below has catalogue e=1 and is absent from the 67. The initial all-catalogue scan selected each for >5 formal σ and >1 km/s. The 120-second replay used identical astrometry and tested the same unconstrained force model. All four returned the correct identity, station/time pairs and inbound barycentric epoch at more than 250 au.

| SOHO object | Excess km/s | (e−1)/σe | Used / returned | Used arc, days | Coordinate RMS, arcsec |
|---|---:|---:|---:|---:|---:|
| C/2004 D4 | 47.326 | 7.56 | 49 / 49 | 1.338 | 98.36 |
| C/2005 D3 | 133.140 | 13.50 | 21 / 21 | 0.662 | 68.18 |
| C/2006 S8 | 72.633 | 13.31 | 34 / 34 | 0.717 | 40.10 |
| C/2006 S9 | 325.345 | 73.87 | 15 / 19 | 0.433 | 43.94 |

These formal excesses coexist with 40–98 arcsec coordinate residuals after fitting. No accepted discovery follows from the unconstrained solutions. [Four follow-up records](soho-index.csv) · [Separate workflow run](https://github.com/dovvnloading/Open-Space/actions/runs/44036442066).

## Inputs and replay

The frozen 2,290-record JPL SBDB query is in [catalog/sbdb.json.gz](../../catalog/sbdb.json.gz). Its uncompressed SHA-256 is `90e391ede544a92c4648436656f38a42cc81f7d2ef7fea3b70dd832c24a26022`. The selection uses `sqrt(GM_sun × max(e−1, 0) / q)` with q in km and the Sun's gravitational parameter in km³/s². The reported barycentric excess speed uses the total solar-system gravitational parameter, approximated by 1.001341838645975 times the solar value. All 62 non-empty bundles came from the same MPC CmtObs gzip, whose source reported September 11, 2026 as its last modification and November 17, 2026 at 06:23 UTC as its retrieval. Source archive SHA-256 is `680c9ee88ef6f615c30ca4faa347a77d5bed5787c8bbc74d8c7257c62a2aa54f`, 25,344,602 bytes. API requests, response hashes or HTTP errors for the remaining five are preserved in the relevant shard records. [Full source query and receipts](evidence/shard-0.json) · [MPC source workflow](https://github.com/dovvnloading/Open-Space/actions/runs/44029645452).

The positional inputs are preserved exactly in `evidence/[supplement-|control-]observations-N.json.gz`, keyed by string JPL SPK ID. Each corresponding `[supplement-|control-]shard-N.json` records hashes of the raw and staged object bytes, the executable and the kernel. All jobs used the same Find_Orb executable SHA-256 `7c126536dafd0c8f42008ba26a6d143624f6626993e0bd523c67603b592699ba` and DE430t SHA-256 `dc8fa8136e71df25c1e4183cd777e6ae536eac099dc6225ad77192cf2192221d`. JPL perihelion time only sets the common output-epoch rule. The six-dimensional fit starts from the astrometry, without a catalogue orbit seed. The `#interstellar` directive removes the high-e solution penalty. The fitted force model includes the planetary perturbers selected by `7fe` and uses gravity alone. [Worker](replay.py) · [Pinned build](../build.sh) · [Six-slice main run](https://github.com/dovvnloading/Open-Space/actions/runs/44035408859).

Run `python verify.py` in this directory to verify all accepted IDs, the source and staged hashes, both observatory cross-checks and the catalogue membership of the two samples using the committed inputs.

To reproduce a whole slice with the source artifact unpacked into `source/`, run from this directory. Python 3.11 or later and the fitter/kernel installed by the pinned build are required.

```sh
python replay.py source 0 output
python replay.py source 0 soho-output --parts 2 --targets soho-targets.txt
python replay.py source 0 control-output --parts 3 --targets iso-targets.txt --join-controls
```

To replay one archived target without the upstream artifact, load the raw lines and fit its frozen row. Set `join_controls=True` only for the separately merged controls.

```python
import gzip, json
from pathlib import Path
from replay import cat, run_one

spkid = "1004106"
rows, _, _ = cat()
row = next(r for r in rows if str(r["spkid"]) == spkid)
bundle = json.loads(gzip.decompress(Path("evidence/observations-5.json.gz").read_bytes()))
result, raw = run_one(row, bundle[spkid].encode("ascii").splitlines(keepends=True))
print(result)
```

## What remains unresolved

The per-object soft ceiling was 120 CPU seconds with a 122-second hard limit and a 300-second wall limit. There was no per-object elapsed-time instrument. The boundary approximation uses a common perihelion-minus-300-years epoch, clipped to stay inside the kernel; a conic at the returned barycentric epoch supplies radius and radial direction. The 100- and 250-au gates have explicit separate columns. A finite formal `sigma_e` and propagated sign under this model cannot substitute for independently validated covariance, nongravitational sensitivity or a reconstructed pre-fragmentation parent orbit. Poor or missing astrometry and the known-control arc losses prevent any population-wide exclusion. The shorter CPU-budget census and the short-arc tracklet injection experiments have their own denominators in the [main study](../../README.md).
