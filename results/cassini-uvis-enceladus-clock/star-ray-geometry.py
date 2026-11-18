import pathlib,sys,math,json,hashlib,numpy as np
p=pathlib.Path(r'D:\Open-Space\campaigns\cassini-uvis-enceladus-20261118');sys.path.insert(0,str(p/'packages'));import spiceypy as s
files=['naif0012.tls','cas00172.tsc','sat427.bsp','200128RU_SCPSE_17034_17098.bsp','pck00010.tpc','enceladus_ssd_250728_v2.tpc']
for f in files:s.furnsh(str(p/f))
ra0=15*(6+58/60+37.54876/3600);dec0=-(28+58/60+19.5102/3600); pmra=3.24;pmdec=1.33
R=252.1
# J2000 catalog position propagated with proper motion to 2017.24; both ray and target use inertial geometric J2000 directions. 
year=2017.235;dec=dec0+(year-2000)*pmdec/3600000;ra=ra0+(year-2000)*pmra/(3600000*math.cos(math.radians(dec0)));u=s.radrec(1,math.radians(ra),math.radians(dec))
def calc(t,abcorr='CN'):
 et=s.str2et(t) if isinstance(t,str) else t; state,lt=s.spkezr('ENCELADUS',et,'J2000',abcorr,'CASSINI');v=state[:3];proj=np.dot(v,u);off=v-proj*u; return {'et':et,'altitude_km':float(np.linalg.norm(off)-R),'distance_km':float(np.linalg.norm(v)),'light_time_s':float(lt),'closest_ray_km':float(proj)}
results=[]
for minute in range(39,46):
 t=f'2017-03-27T14:{minute:02}:00';g=calc(t); h=calc(g['et']+1);g['utc']=t;g['one_second_altitude_change_km']=h['altitude_km']-g['altitude_km'];g.pop('et');results.append(g)
out={'star':{'name':'epsilon CMa HIP 33579','source':'https://simbad.u-strasbg.fr/simbad/sim-id?Ident=eps+Cma&submit=submit+id','ICRS_J2000_ra_deg':ra0,'ICRS_J2000_dec_deg':dec0,'pmra_cosdec_mas_year':pmra,'pmdec_mas_year':pmdec,'propagated_epoch':year,'propagated_ra_deg':ra,'propagated_dec_deg':dec},'model':'J2000 geometric infinite-distance stellar ray with proper motion; apparent Enceladus center by converged Newtonian light-time CN, relative to Cassini reception time; no stellar aberration because both directions are in inertial geometric frame; 252.1 km spherical radius; no terrain or atmosphere. Stellar parallax, aberration-frame residual, and star parameter uncertainty unquantified.','kernels':[{'filename':f,'sha256':hashlib.sha256((p/f).read_bytes()).hexdigest()} for f in files],'samples':results}
(p/'star-ray-geometry.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
