"""Extract the previously published row for the distance-sensitive BOSS target."""
import csv
from hashlib import sha256
from urllib.request import urlopen
from astropy.io import fits
from fetch_survey import HERE,WORK
URL='https://data.sdss.org/sas/dr19/vac/mwm/minesweeper/minesweeper_v1.0.0.fits'
SHA='8afae40c339655bb0e8c4f37ac9faddd460d0aaa0973dfb82326313becad314b'
if __name__=='__main__':
    path=WORK/'minesweeper_v1.0.0.fits'
    if not path.is_file():
        with urlopen(URL,timeout=180) as response:path.write_bytes(response.read())
    assert sha256(path.read_bytes()).hexdigest()==SHA
    cols=('source_id sdssid snr n_spall mjd Dist Dist_lerr Dist_uerr Dist_err Vrad Vrad_err '
          'Teff logg logg_err FeH Av Mass R_gal Vx_gal Vy_gal Vz_gal E_kin_mw22 E_pot_mw22 '
          'E_tot_mw22 E_tot_mw22_err R_peri_mw22 FLAG in_substructure in_cluster').split()
    output=HERE/'audit'/'minesweeper.csv';output.parent.mkdir(exist_ok=True)
    with fits.open(path) as hdus, output.open('w',newline='') as fp:
        a=hdus[1].data
        matched=a[a['source_id']==4406622293782568192]
        if len(matched)!=1:raise ValueError('MINESweeper row missing or ambiguous')
        writer=csv.writer(fp);writer.writerow(cols);writer.writerow([matched[0][c] for c in cols])
