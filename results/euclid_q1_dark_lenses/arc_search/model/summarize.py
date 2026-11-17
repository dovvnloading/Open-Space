"""Collect five spatial folds, including optional broad-centre fits, into reproducible tables."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def write(path,rows):
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=rows[0].keys(),lineterminator="\n")
        w.writeheader();w.writerows(rows)


def summarize(folds,output):
    entries={}
    for directory in [output,folds]:
        for path in directory.glob("*.json"):
            r=json.loads(path.read_text())
            if "heldout_pixels" in r:
                entries[(r["id"],r.get("center_bound_arcsec",.35),r.get("fold",0))]=r
    data=[]
    summary=[]
    for name,bound in sorted(set((key[0],key[1]) for key in entries)):
        fits=[entries[(name,bound,i)] for i in range(5)]
        for i,r in enumerate(fits):
            l,k=r["lens"],r["knots"]
            data.append({"id":name,"center_bound_arcsec":bound,"fold":i,
                         "heldout_n":r["heldout_pixels"],"training_n":r["training_pixels"],
                         "test_lens":l["test_mean_sq"],"test_knots":k["test_mean_sq"],
                         "train_lens":l["train_mean_sq"],"train_knots":k["train_mean_sq"],
                         "b_arcsec":l["parameters"][0],"gamma1":l["parameters"][3],"gamma2":l["parameters"][4],
                         "lens_dx_arcsec":l["parameters"][1],"lens_dy_arcsec":l["parameters"][2],
                         "lens_parameters":json.dumps(l["parameters"]),"knot_parameters":json.dumps(k["parameters"]),
                         "lens_amplitudes_native":json.dumps(l["amplitudes_native"]),
                         "knot_amplitudes_native":json.dumps(k["amplitudes_native"])})
        n=np.array([r["heldout_pixels"] for r in fits])
        lens=np.average([r["lens"]["test_mean_sq"] for r in fits],weights=n)
        knots=np.average([r["knots"]["test_mean_sq"] for r in fits],weights=n)
        summary.append({"id":name,"center_bound_arcsec":bound,"ra_deg":fits[0]["ra"],"dec_deg":fits[0]["dec"],
                        "heldout_pixels_total":int(n.sum()),"radius_proposal":fits[0]["proposal_radius_arcsec"],
                        "number_knots":len(fits[0]["detected_knots"]),"mean_square_lens":f"{lens:.5f}",
                        "mean_square_knots":f"{knots:.5f}","delta_lens_minus_knots":f"{lens-knots:.5f}",
                        "folds_favour_knots":sum(r["test_delta_lens_minus_knots"]>0 for r in fits),
                        "b_range_arcsec":json.dumps([round(float(v),4) for v in
                                                     np.percentile([r["lens"]["parameters"][0] for r in fits],[0,100])]),
                        "folds_center_hits_bound":sum(max(abs(v) for v in r["lens"]["parameters"][1:3]) >= bound-.005
                                                       for r in fits)})
    write(output/"cross_validation_folds.csv",data)
    write(output/"cross_validation.csv",summary)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folds",type=Path)
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    args=p.parse_args()
    summarize(args.folds,args.output)
