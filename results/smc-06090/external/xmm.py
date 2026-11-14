"""Check target and matched control sources in twelve XMM-Newton OM source lists."""
import csv,json,hashlib
from astropy.io import fits
from astropy.time import Time
from scipy.spatial import cKDTree
from common import D,np,SKY,hjd,shape
def xyz(r):
    return np.column_stack(((r['RA_CORR']-SKY.ra.deg)*np.cos(SKY.dec.radian),r['DEC_CORR']-SKY.dec.deg))*3600
all_rows=[];files=[];m2=[]
for path in sorted((D/'inputs').glob('P*FTZ')):
    with fits.open(path) as hd:
        h=hd[0].header;r=hd[1].data.copy()
    pos=xyz(r); ix=np.argmin(np.linalg.norm(pos,axis=1))
    ra=r['RA_CORR'][ix];dec=r['DEC_CORR'][ix]
    t0,t1=Time([h['DATE-OBS'],h['DATE-END']],scale='utc').mjd;j=hjd(.5*(t0+t1))
    ph,d1,d2=shape(j);r0=r[ix]
    all_rows.append(dict(file=path.name,filter=h['FILTER'],start=h['DATE-OBS'],end=h['DATE-END'],hjd_mid=j,exposure_s=h['EXPOSURE'],phase=ph,delta_primary_cycles=d1,delta_secondary_cycles=d2,ra=float(ra),dec=float(dec),sep_arcsec=float(np.linalg.norm(pos[ix])),mag=float(r0['MAG']),err=float(r0['MAGERR']),rate=float(r0['CORR_RATE']),rateerr=float(r0['CORR_RATE_ERR']),qflag=int(r0['QFLAG']),cflag=int(r0['CFLAG'])))
    files.append(dict(name=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),url='https://nasa-heasarc.s3.amazonaws.com/xmm/data/rev0/0741450101/PPS/'+path.name))
    if h['FILTER']=='UVM2':m2.append((r,pos,ix))
with (D/'xmm-photometry.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=all_rows[0]);w.writeheader();w.writerows(all_rows)
rr,xy,k=m2[1]
controls=[]
for m,row in enumerate(rr):
    if k==m or not 15.2<row['MAG']<17.7 or np.linalg.norm(xy[m])>360:continue
    settings=[]
    for r,pos,_ in m2:
        dist,i=cKDTree(pos).query(xy[m])
        if dist>1 or r['QFLAG'][i]!=0 or r['MAGERR'][i]>.14:break
        settings.append([float(r['MAG'][i]),float(r['MAGERR'][i])])
    if len(settings)==3:
        mag,err=np.array(settings).T
        shift=mag[1]-.5*(mag[0]+mag[2])
        controls.append([m,float(shift),float(np.sqrt(err[1]**2+.25*(err[0]**2+err[2]**2)))])
md=np.median([x[1] for x in controls]);mad=1.4826*np.median(np.abs(np.array([x[1] for x in controls])-md))
r=[x for x in all_rows if x['filter']=='UVM2'];delta=r[1]['mag']-.5*(r[0]['mag']+r[2]['mag'])
out=dict(files=files,obsid='0741450101',primary_midpoint_offset_min=r[1]['delta_primary_cycles']*1.1696719*1440,delta_mag=delta,delta_mag_stat=np.sqrt(r[1]['err']**2+.25*(r[0]['err']**2+r[2]['err']**2)),n_controls=len(controls),control_median_mag=md,control_robust_scatter_mag=mad,control_tail=sum(x[1]>=delta for x in controls),controls=controls)
(D/'xmm.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k not in ['files','controls']})
