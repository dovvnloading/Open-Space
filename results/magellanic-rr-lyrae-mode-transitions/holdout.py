"""Compare OGLE-only clock predictions on independent data and stress the blend test."""
import argparse
import json
from pathlib import Path
import numpy as np

from checks import flux_check
from measure import NAME, SOURCE, gaia
from screen import model, photo


def ratio_predict(d, periods):
    """Leave every visit of one night out, with an independently supplied fixed period."""
    t, y, e = d.T
    X = model(t, periods, 2)
    w = 1/np.maximum(.015,e)
    predictions = np.zeros(len(t))
    for night in np.unique(np.floor(t)):
        test = np.floor(t)==night
        keep = ~test
        b = np.linalg.lstsq(X[keep]*w[keep,None], y[keep]*w[keep], rcond=None)[0]
        predictions[test] = X[test] @ b
    return float(np.sqrt(np.mean((y-predictions)**2)))


def main(cache):
    rng=np.random.default_rng(902211)
    q=f"SELECT * FROM gaiadr3.epoch_photometry WHERE source_id={SOURCE}"
    row=gaia(q,cache)[0][0]
    def col(k):
        return np.array(row[k].split(),dtype=float)
    d=np.c_[col("g_transit_time")+5197.5,col("g_transit_mag"),
            (2.5/np.log(10))/col("g_transit_flux_over_error")]
    d=d[(np.array(row["variability_flag_g_reject"].split())=="False")&np.isfinite(d).all(axis=1)]
    src=json.loads((Path(__file__).parent/"measurements.json").read_text())["windows"]
    early=src["4I_5750_6120"]["periods"]
    late=src["4I_6850_7215"]["periods"]
    tested={"constant":[], "early_overtone":[early[0]],"early_fundamental":[early[1]],
            "late_fundamental":[late[1]]}
    gaia_rms={k:ratio_predict(d,periods) for k,periods in tested.items()}
    # Independent late OGLE-I season. All tested periods come from earlier OGLE-I seasons.
    iv=photo(4,NAME,"I",cache)
    last=iv[(iv[:,0]>=7215)&(iv[:,0]<7580)]
    last_rms={k:ratio_predict(last,periods) for k,periods in tested.items()}
    # Allow the per-mode waveform to vary across halves of the best populated season.
    mid=iv[(iv[:,0]>=5750)&(iv[:,0]<6120)]
    halves=np.array_split(mid[mid[:,0].argsort()],2)
    mixed=[flux_check(part,early,rng,5000) for part in halves]
    return dict(ogle_periods=tested,gaia_n=len(d),gaia_night_loo_rmse=gaia_rms,
                independent_late_I_n=len(last),independent_late_I_night_loo_rmse=last_rms,
                middle_season_half_flux_tests=mixed)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=Path("raw-cache"))
    parser.add_argument("--out",type=Path,default=Path("holdout.json"))
    args=parser.parse_args()
    args.out.write_text(json.dumps(main(args.cache),indent=2)+"\n")
