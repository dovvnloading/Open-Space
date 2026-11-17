"""Within-visit scatter of the Hephaistos A–J NEOWISE final-year detections."""
import argparse
import csv
import math
from pathlib import Path
from scipy.special import gammaincc
from measure import single_exposures
from fetch import targets


def compute(archive):
    records = []
    for target in targets():
        for band in (1,2):
            visits = []
            for point in single_exposures(target, band, archive):
                if not visits or point[1] - visits[-1][-1][1] > 20:
                    visits.append([])
                visits[-1].append(point)
            for visit in visits:
                if len(visit) < 3:
                    continue
                # Independent-point test is descriptive for same-scan exposures.
                weights = [1/(point[3]**2 + .02**2) for point in visit]
                mag = sum(w*p[2] for w,p in zip(weights,visit)) / sum(weights)
                chi2 = sum(w*(p[2]-mag)**2 for w,p in zip(weights,visit))
                dof = len(visit)-1
                records.append(dict(label=target['label'],band=f'W{band}',mjd_start=visit[0][1],
                                    n=len(visit),weighted_mag=mag,chi2=chi2,dof=dof,
                                    p_raw=gammaincc(dof/2,chi2/2)))
    for rec in records:
        rec['p_bonferroni'] = min(1.,len(records)*rec['p_raw'])
    return records


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('archive',type=Path,help='Extracted archive_inputs.zip')
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    records=compute(args.archive)
    with args.output.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,records[0].keys(),lineterminator='\n')
        writer.writeheader()
        for r in records:
            for name in ('mjd_start','weighted_mag','chi2','p_raw','p_bonferroni'):
                r[name]=f'{r[name]:.8g}' if name!='mjd_start' else f'{r[name]:.8f}'
            writer.writerow(r)
