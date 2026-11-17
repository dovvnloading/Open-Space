"""Fetch the Gaia DR3 rows by integer ID and preserve the exact received CSV bytes."""
import argparse
import csv
import hashlib
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ENDPOINT='https://gea.esac.esa.int/tap-server/tap/sync'
COLS=('source_id ra dec parallax pmra pmdec ra_error dec_error parallax_error pmra_error pmdec_error '
      'ra_dec_corr ra_parallax_corr ra_pmra_corr ra_pmdec_corr dec_parallax_corr dec_pmra_corr dec_pmdec_corr '
      'parallax_pmra_corr parallax_pmdec_corr pmra_pmdec_corr ruwe visibility_periods_used duplicated_source '
      'astrometric_params_solved nu_eff_used_in_astrometry pseudocolour ecl_lat phot_g_mean_mag bp_rp '
      'radial_velocity radial_velocity_error rv_nb_transits rv_chisq_pvalue rv_renormalised_gof non_single_star '
      'ipd_frac_multi_peak ipd_gof_harmonic_amplitude').split()


def fetch(args):
    i, ids, out = args
    sql = 'SELECT '+','.join(COLS)+' FROM gaiadr3.gaia_source WHERE source_id IN ('+','.join(ids)+') ORDER BY source_id'
    qp=urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql))
    path=out/f'rows_{i:03}.csv'
    if not path.exists():
        for attempt in range(4):
            try:
                raw=urlopen(ENDPOINT+'?'+qp,timeout=120).read()
                table=list(csv.DictReader(io.StringIO(raw.decode())))
                if not table or not {r['source_id'] for r in table}<=set(ids):
                    raise ValueError('empty or invalid Gaia result')
                path.write_bytes(raw)
                break
            except Exception:
                if attempt==3: raise
                time.sleep(2**attempt)
    rows=list(csv.DictReader(path.open()))
    return dict(query=sql,rows=len(rows),missing_ids=sorted(set(ids)-{r['source_id'] for r in rows}),
                file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('out',type=Path);p.add_argument('--key',default='SOURCE_ID');a=p.parse_args()
    ids=sorted({str(int(r[a.key])) for r in csv.DictReader(open(a.input))},key=int)
    a.out.mkdir(exist_ok=True,parents=True)
    with ThreadPoolExecutor(max_workers=4) as ex:
        receipt=list(ex.map(fetch,[(i//60,ids[i:i+60],a.out) for i in range(0,len(ids),60)]))
    (a.out/'receipt.json').write_text(json.dumps({'endpoint':ENDPOINT,'blocks':receipt},indent=2)+'\n')
