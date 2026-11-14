"""Compare the orbital period with a published sample and compute a conditional GR scale."""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import subprocess
from urllib.request import urlretrieve

ROOT=Path(__file__).resolve().parent
COLUMNS=('ogle_iii_id','sidereal_period_days','periastron_period_years')
def parse_paper(text):
    table=text.split('Apsidal motion elements of 90 EBs')[-1].split('Notes−Listed parameters')[0]
    rows=[]
    for line in table.splitlines():
        cells=line.split()
        if not cells or not re.fullmatch(r'SMC-ECL-\d{4}',cells[0]):
            continue
        if len(cells)<9:raise ValueError(line)
        rows.append((cells[0],float(cells[2].split('(')[0]),float(cells[-1].split('(')[0])))
    if len(rows)!=90:raise ValueError(f'Expected 90 published targets, got {len(rows)}')
    return rows

def main():
    cmd=argparse.ArgumentParser(description=__doc__)
    cmd.add_argument('--table-text',type=Path)
    cmd.add_argument('--download',action='store_true')
    opt=cmd.parse_args()
    if opt.download:
        path=ROOT/'hong2016.pdf';urlretrieve('https://arxiv.org/pdf/1604.02200',path)
        opt.table_text=ROOT/'hong2016.txt'
        subprocess.run(['pdftotext','-layout',path,opt.table_text],check=True)
    if opt.table_text:
        rows=parse_paper(opt.table_text.read_text())
        with (ROOT/'hong2016_table4.csv').open('w') as f:
            w=csv.writer(f);w.writerow(COLUMNS);w.writerows(rows)
    with (ROOT/'hong2016_table4.csv').open() as f:
        comparison=list(csv.DictReader(f))
    joint=json.loads((ROOT.parent/'joint'/'fits.json').read_text())['full']['kepler']
    period=1.1696719
    u=joint['period_years']; e=joint['edge_on_eccentricity']
    observed=360*period/(365.25*u)  # degrees per sidereal revolution
    # Gimenez (1985); mass in solar units, period in days, rate in deg/cycle.
    def gr(mass):return 5.447e-4/(1-e**2)*(mass/period)**(2/3)
    req=period*(observed*(1-e**2)/5.447e-4)**1.5
    result=dict(target_p_days=period,target_ogle_iii_I=18.232,hong2016_magnitude_cut='I < 18.0',
                table_rows=len(comparison),shorter=[r for r in comparison if float(r[COLUMNS[1]])<period],
                joint_model_cycle_years=u,conditional_e=e,forward_rate_deg_per_orbit=observed,
                total_mass_solar_for_gr_alone=req,
                gr_examples={str(m):{'period_years':360*period/(365.25*gr(m)),'fraction_of_candidate_rate':gr(m)/observed}
                             for m in [2,5,10,20,50]})
    (ROOT/'context.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)
if __name__=='__main__':main()
