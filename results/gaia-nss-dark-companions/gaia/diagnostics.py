"""Retrieve Gaia DR3 variability classes, blend diagnostics and Gaia RV summaries for the screen."""
import argparse
import csv
import hashlib
import io
import json
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ENDPOINT = 'https://gea.esac.esa.int/tap-server/tap/sync'
GAIA = ('source_id ra dec ref_epoch pmra pmdec parallax parallax_error ruwe '
        'ipd_frac_multi_peak ipd_gof_harmonic_amplitude ipd_frac_odd_win '
        'visibility_periods_used duplicated_source phot_g_mean_mag bp_rp phot_bp_rp_excess_factor '
        'phot_variable_flag radial_velocity radial_velocity_error rv_nb_transits rv_nb_deblended_transits '
        'rv_time_duration rv_amplitude_robust rv_chisq_pvalue rv_renormalised_gof '
        'rv_expected_sig_to_noise has_rvs has_xp_continuous').split()
VAR = ['source_id','best_class_name','best_class_score']
AP = ('source_id teff_gspphot logg_gspphot mh_gspphot ag_gspphot logposterior_msc flags_msc '
      'teff_msc1 teff_msc2 mass_flame mass_flame_lower mass_flame_upper flags_flame').split()
TABLES = [('gaia_source',GAIA),('vari_classifier_result',VAR),('astrophysical_parameters',AP)]


def query(task):
    table, fields, batch, cache = task
    adql = ('SELECT '+','.join(fields)+' FROM gaiadr3.'+table+' WHERE source_id IN ('+
            ','.join(batch)+') ORDER BY source_id')
    digest = hashlib.sha256(adql.encode()).hexdigest()
    dest = cache / (digest + '.csv')
    if not dest.exists():
        url = ENDPOINT + '?' + urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',
                                                          MAXREC=len(batch)+1,QUERY=adql))
        data = urllib.request.urlopen(url,timeout=120).read()
        if data.splitlines()[0].decode().casefold() != ','.join(fields):
            raise ValueError(data[:200])
        dest.write_bytes(data)
    data = dest.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    ids = {r['source_id'] for r in rows}
    if len(rows) != len(ids) or not ids.issubset(batch) or (table=='gaia_source' and ids!=set(batch)):
        raise ValueError('missing or unexpected identifiers for ' + table)
    return rows,dict(table=table,query_sha256=digest,response_sha256=hashlib.sha256(data).hexdigest(),
                     requested=len(batch),received=len(rows))


def run(directory, cache):
    cache.mkdir(parents=True,exist_ok=True)
    with (directory/'mass_screen.csv').open() as file:
        ids = sorted({r['source_id'] for r in csv.DictReader(file)},key=int)
    records = {i:{'source_id':i} for i in ids}
    tasks = [(table,fields,ids[i:i+120],cache) for table,fields in TABLES for i in range(0,len(ids),120)]
    manifest = []
    with ThreadPoolExecutor(3) as pool:
        for rows,m in pool.map(query,tasks):
            manifest.append(m)
            for r in rows:
                records[r['source_id']].update(r)
    with (directory/'source_diagnostics.csv').open('w', newline='') as f:
        w=csv.DictWriter(f,GAIA+VAR[1:]+AP[1:],lineterminator='\n')
        w.writeheader()
        w.writerows(records.values())
    (directory/'diagnostic_inputs.json').write_text(json.dumps(dict(endpoint=ENDPOINT,chunks=manifest),
                                                              indent=2)+'\n')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dir', default=Path(__file__).parent,type=Path)
    parser.add_argument('--cache',required=True,type=Path)
    args=parser.parse_args()
    run(args.dir,args.cache)
