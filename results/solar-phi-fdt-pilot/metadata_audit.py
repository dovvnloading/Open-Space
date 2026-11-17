"""Audit scan timing and processing-table identifiers without image arrays."""
import argparse
import json
from pathlib import Path
from astropy.io import fits
from astropy.time import Time
from geometry import read_header


def audit(path):
    h=read_header(path)
    # These two small binary-table extensions contain acquisition metadata,
    # not the primary image array. Do not access HDU 0 .data.
    with fits.open(path,memmap=True,do_not_scale_image_data=True) as hdus:
        exp=hdus['PHI_PERIMAGE_information'].data
        stamps=[str(x).strip() for x in exp['EXP_START_TIME']]
        if len(stamps)!=len(set(stamps)):raise ValueError('duplicate image timestamps')
        times=Time(stamps,format='isot',scale='utc')
        if any(t<Time(h['DATE-BEG']) or t>Time(h['DATE-END']) for t in times):raise ValueError('timestamp outside scan')
        blk=hdus['PHI_FITS_PROC_block'].data
        ids=sorted({int(x) for x in blk['PHI_PROC_operandDID'] if int(x)!=0})
        accum=hdus['PHI_FITS_ACCU_settings'].data
        counts=sorted({int(x) for x in accum['PHI_ACCU_accumulations']})
    return {'filename':h['FILENAME'],'date_beg':h['DATE-BEG'],'date_end':h['DATE-END'],'perimage_records':len(stamps),'first':min(stamps),'last':max(stamps),'accumulations':counts,'operand_ids':ids,'timestamps':stamps}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('files',nargs='+',type=Path);a=parser.parse_args()
    entries=[audit(p) for p in a.files]
    for i in range(len(entries)):
        for j in range(i):
            if set(entries[i]['timestamps'])&set(entries[j]['timestamps']):raise ValueError('shared per-image timestamp')
            if max(Time(entries[i]['date_beg']),Time(entries[j]['date_beg']))<min(Time(entries[i]['date_end']),Time(entries[j]['date_end'])):raise ValueError('overlapping scan windows')
    for e in entries:e.pop('timestamps')
    print(json.dumps(entries,indent=2))
