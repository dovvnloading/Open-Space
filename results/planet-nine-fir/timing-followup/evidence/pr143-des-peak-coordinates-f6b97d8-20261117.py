import json,csv,gzip,numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
base=Path('D:/Open-Space/campaigns');hdr=fits.getheader(base/'pr143-des-inputs-581c9c4-20261117/DECam_605285.response');w=WCS(hdr)
for name in ('dec1-f6b97d8-296a210','dec2-f6b97d8-b000d53'):
 p=base/('pr143-des-'+name+'-20261117')/'run'; peak=json.loads((p/'dense_peaks.json').read_text())[0];track=list(csv.DictReader(gzip.open(p/'difference_tracks.csv.gz','rt')));x,y=w.world_to_pixel_values(peak['ra'],peak['dec']);nearest=track[peak['nearest_track']];tx,ty=w.world_to_pixel_values(float(nearest['ra_a']),float(nearest['dec_a']));print(name,'feature native WCS pixel x,y',x,y,'nearest',peak['nearest_track'],'track native',tx,ty,'distance',np.hypot(x-tx,y-ty),'peak vector',peak['dx'],peak['dy'],'nearest valid',nearest['valid'],'track flag',nearest['flagged'])
