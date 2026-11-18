"""Hash mission/release descriptions that scope external time-domain coverage."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from fetch import OUT, get, write_csv, SHA
URLS={
'UKIRT WFIRST precursor 2015-19':'https://exoplanetarchive.ipac.caltech.edu/docs/UKIRTMission.html',
'HST THOR GO-17776 early release':'https://archive.stsci.edu/hlsp/thor',
'HST GO-17776 scientific description':'https://arxiv.org/abs/2605.06778',
'Euclid Q2 release, IRSA':'https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_q2.html',
'Euclid Q2 temporal sampling':'https://arxiv.org/abs/2606.25883',
'Roman planned observations':'https://roman-docs.stsci.edu/roman-community-defined-surveys/roman-observations-in-the-first-two-years-of-science-operations',
'ESA Gaia end of observations':'https://www.esa.int/Science_Exploration/Space_Science/Gaia/Last_starlight_for_ground-breaking_Gaia',
'Gaia public light curve definitions':'https://gsaweb.ast.cam.ac.uk/alerts/docs/',
'GaiaX per-run detections':'https://gsaweb.ast.cam.ac.uk/alerts/gaiax/',
'PRIME bulge survey description':'https://arxiv.org/abs/2508.14474',
'Historical MOA-II statistical archive':'https://exoplanetarchive.ipac.caltech.edu/docs/MOAMission.html',
'ROME/REA archival data release':'https://exoplanetarchive.ipac.caltech.edu/docs/rome/ROME_home.html'}
def run(kv):
    name,url=kv
    try:
        data=get(url);return dict(source=name,url=url,sha256=SHA(data),bytes=len(data),status='ok')
    except Exception as e:return dict(source=name,url=url,status='error',error=str(e))
if __name__=='__main__':
    with ThreadPoolExecutor(6) as ex:write_csv(OUT/'archive_receipts.csv',list(ex.map(run,URLS.items())))
