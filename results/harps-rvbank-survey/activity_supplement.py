"""Match published Mount Wilson S-indices to individual southern HARPS exposures.

Ibáñez Bustos et al. (2025), DOI 10.1051/0004-6361/202450348,
CDS J/A+A/696/A230/tableb1. The catalog contains one instrument and MJD
per spectrum. Compare this time to the HARPS UTC filename timestamp.
"""
import json
import re

import numpy as np
import pandas as pd

TABLE = "J/A+A/696/A230/tableb1"
FIELDS = {"name": ("Name", 1, 9), "mjd": ("MJD", 11, 26),
          "inst": ("Inst", 30, 37), "program": ("ProgID", 39, 52),
          "s": ("S-Index", 56, 67), "comment": ("Com", 69, 98)}
POSITIONS = {"name": ("Name", 1, 9), "rah": ("RAh", 12, 13),
             "ram": ("RAm", 15, 16), "ras": ("RAs", 18, 22),
             "sign": ("DE-", 24, 24), "ded": ("DEd", 25, 26),
             "dem": ("DEm", 28, 29), "des": ("DEs", 31, 34)}


def star_id(value):
    s = re.sub(r"\s+", "", str(value)).upper()
    # Astrometry and time checks are independent of this spelling normalization.
    return re.sub(r"^(?:GL|GJ|GLIESE)[_-]?0*(?=\d)", "GJ", s)


def add_published_sindex(data, work, out, manifest, errors, catalog):
    """Use one S-index reduction per star. Leave ambiguous matches unassigned."""
    details = {"catalog": TABLE, "status": "fetching"}
    matched = []
    fields = ["star", "source_name", "epoch_id", "source_mjd", "source_program",
              "bjd", "sindex", "utc_difference_seconds", "identity_method"]
    try:
        aux = catalog(TABLE, FIELDS, 2808, work, manifest, errors)
        aux["name"] = aux.name.map(star_id)
        aux["inst"] = aux.inst.astype("string").str.strip().str.upper()
        aux["s"] = pd.to_numeric(aux.s, errors="coerce")
        aux["mjd"] = pd.to_numeric(aux.mjd, errors="coerce")
        details["all_instruments"] = aux.inst.value_counts().to_dict()
        # The published comments can describe flares and other exclusions.
        clean = aux.comment.fillna("").astype(str).str.strip().isin(["", "-", "--"])
        aux = aux[aux.inst.eq("HARPS") & aux.mjd.between(52700, 65000) &
                  np.isfinite(aux.s) & aux.s.gt(0) & clean].copy()
        details["harps_clean_positive_s"] = len(aux)
        used = data.loc[data.disposition.eq("used")].copy()
        saved = set(used.loc[used.sindex.notna(), "star"])
        eligible = used[~used.star.isin(saved)].copy()
        eligible["name"] = eligible.star.map(star_id)
        # A UTC MJD cannot be compared directly with barycentric BJD.
        stamps = eligible.epoch_id.astype("string").str.extract(
            r"(\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}\.\d{3})(?!\d)", expand=False)
        utc = pd.to_datetime(stamps, utc=True, errors="coerce", format="ISO8601")
        eligible["utc_mjd"] = (utc - pd.Timestamp("1858-11-17", tz="UTC")).dt.total_seconds()/86400
        eligible = eligible[eligible.utc_mjd.notna()]
        positions = None
        unmatched = set(aux.name) - set(eligible.name)
        if unmatched:
            try:
                pos = catalog("J/A+A/696/A230/sources", POSITIONS, 35, work, manifest, errors)
                decsign = pos.sign.astype("string").str.strip()
                for field in ("rah","ram","ras","ded","dem","des"):
                    pos[field] = pd.to_numeric(pos[field], errors="coerce")
                pos["ra"] = 15 * (pos.rah + pos.ram/60 + pos.ras/3600)
                pos["dec"] = (pos.ded + pos.dem/60 + pos.des/3600) * decsign.map({"-": -1, "+": 1})
                pos["name"] = pos.name.map(star_id)
                positions = pos[~pos.name.duplicated(keep=False)].set_index("name")
            except Exception as exc:
                errors.append(f"{TABLE}: position fallback unavailable: {exc}")
        # Include already-covered targets when testing uniqueness on the sky.
        sky = used.groupby("star")[["ra_deg","dec_deg"]].median()
        ids = used.groupby("star").size().index.to_series().map(star_id)
        associations = {}
        for name in sorted(set(aux.name)):
            hits = list(ids[ids == name].index)
            method = "catalog identifier"
            if not hits and positions is not None and name in positions.index:
                rec = positions.loc[name]
                dra = (sky.ra_deg - rec.ra + 180) % 360 - 180
                sep = np.hypot(dra*np.cos(np.deg2rad(rec.dec)), sky.dec_deg-rec.dec)
                hits = list(sky.index[sep <= 2/3600])
                method = "unique J2000 position within 2 arcsec"
            if len(hits) == 1 and hits[0] not in saved:
                associations[name] = (hits[0], method)
        mapped = pd.Series([v[0] for v in associations.values()], dtype=str)
        ambiguous = set(mapped[mapped.duplicated(keep=False)])
        for name, (star, method) in associations.items():
            if star in ambiguous:
                continue
            left = eligible[eligible.star.eq(star)].sort_values("utc_mjd")
            right = aux[aux.name.eq(name)].sort_values("mjd")
            lt, rt = left.utc_mjd.to_numpy(), right.mjd.to_numpy()
            if not len(lt) or not len(rt):
                continue
            tol = 1.0 / 86400
            l = np.searchsorted(rt, lt-tol)
            h = np.searchsorted(rt, lt+tol, side="right")
            rl = np.searchsorted(lt, rt-tol)
            rh = np.searchsorted(lt, rt+tol, side="right")
            mutual = (h-l==1) & ((rh-rl)[l.clip(0,len(rt)-1)]==1)
            for i in np.flatnonzero(mutual):
                j = l[i]
                a, b = left.iloc[i], right.iloc[j]
                matched.append(dict(star=star, source_name=name, epoch_id=a.epoch_id,
                                    source_mjd=b.mjd, source_program=b.program, bjd=a.bjd,
                                    sindex=float(b.s), utc_difference_seconds=(lt[i]-rt[j])*86400,
                                    identity_method=method))
        details.update(status="completed", stars_with_accepted_matches=len({r["star"] for r in matched}),
                       accepted_matches=len(matched))
        ledger = pd.DataFrame(matched, columns=fields)
        # Only commit new measurements after all identities and joins have succeeded.
        if not ledger.empty:
            chosen = eligible.reset_index(names="row").merge(
                ledger, on=["star", "epoch_id", "bjd"], validate="one_to_one", suffixes=("", "_new"))
            data.loc[chosen.row, "sindex"] = chosen.sindex_new.to_numpy()
    except Exception as exc:
        details.update(status="failed", error=str(exc))
        raise
    finally:
        pd.DataFrame(matched, columns=fields).to_csv(out / "ibanez-sindex-matches.csv", index=False)
        (out / "ibanez-sindex-coverage.json").write_text(json.dumps(details, indent=2)+"\n")
