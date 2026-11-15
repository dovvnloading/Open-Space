"""Tabulate the official YMW16 v1.3 distance-to-DM model at the timed positions.

Pass the extracted ymw16_v1.3 source directory as the sole argument. Compile the
authors' model with a tiny main function to avoid legacy CLI parsing and enable
an ordinary zero-on-success process return. The numerical model is unmodified.
"""
import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from analyze import ROOT, read, sexagesimal

DRIVER = r"""
#include <stdlib.h>
#include "cn.h"
int main(int argc, char **argv) {
    if (argc != 4) return 1;
    char directory[] = "./", label[] = "";
    dmdtau(atof(argv[1]), atof(argv[2]), atof(argv[3]), 100, 2, 1, 0, directory, label);
    return 0;
}
"""
MATRIX = np.array([[-.0548755604, -.8734370902, -.4838350155],
                   [ .4941094279, -.4448296300,  .7469822445],
                   [-.8676661490, -.1980763734,  .4559837762]])


def galactic(row):
    a, d = np.radians([15 * sexagesimal(row['ra_j2000']), sexagesimal(row['dec_j2000'])])
    x, y, z = MATRIX @ [np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)]
    return np.degrees(np.arctan2(y, x)) % 360, np.degrees(np.arcsin(z))


def generate(source):
    modules = ('dmdtau dora fermibubble frb_d galcen gum lmc localbubble ne_crd '
               'nps smc spiral thick thin ymw16par').split()
    with tempfile.TemporaryDirectory() as temp:
        main, executable = Path(temp)/'main.c', Path(temp)/'ymw16'
        main.write_text(DRIVER)
        command = ['cc', '-O2', '-fcommon', '-I', str(source), str(main),
                   *(str(source/(s+'.c')) for s in modules), '-lm', '-o', str(executable)]
        subprocess.run(command, check=True)
        distances = (8700, 9260, 9820, 10600)
        with (ROOT/'smooth-galaxy.csv').open('w') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(['pulsar', 'l_deg', 'b_deg', *(f'dm_{d}pc' for d in distances)])
            for row in read('fast.csv'):
                l, b = (f'{v:.9f}' for v in galactic(row))
                dms = []
                for d in distances:
                    out = subprocess.check_output([str(executable), l, b, str(d)],
                                                  cwd=source, text=True)
                    dms.append(re.search(r'DM:\s*([\d.]+)', out)[1])
                writer.writerow([row['pulsar'], l, b, *dms])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    generate(parser.parse_args().source.resolve())
