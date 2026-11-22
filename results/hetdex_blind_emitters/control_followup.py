"""Run the same spectral vetting and coordinate-based follow-up on shuffled leads.

Each null spectrum retains the assigned source row's metadata and sky position.
Independent visit measurements and ancillary checks use the actual sky at that
position, measuring the possibility of accidental corroboration in this design.
"""
from pathlib import Path

from scan import BATCH, jsonlines
from spectral import detail, shuffle_fibres
from unknown_lines import unknown_spectrum
from vetting import vet


def control_followup(cat,out):
    from itertools import groupby
    root=Path(out)/'shuffle-vet'
    root.mkdir(exist_ok=True)
    def records():
        rows=jsonlines(Path(out)/'shuffled.jsonl.gz')
        for start, subgroup in groupby(rows,lambda r:r['row_index']//BATCH*BATCH):
            f,e=cat.read(start,min(start+BATCH,cat.size))
            fy,fe=shuffle_fibres(f,e,seed=410087+start)
            for row in subgroup:
                i=row['row_index']
                m=cat.metadata(i)
                fits=detail(cat.wave,fy[i-start],fe[i-start],row,density_unit=1e-17)
                if row['channel']=='unknown':
                    fits['unknown_spectrum']=unknown_spectrum(cat.wave,fy[i-start],fe[i-start],row,density_unit=1e-17)
                yield {**m, 'catalog':m, **row, **fits, 'shuffled_detail_complete':True}
    vet(records(),cat,root)
