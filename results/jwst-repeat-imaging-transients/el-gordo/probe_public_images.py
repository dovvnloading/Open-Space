import csv,io
from pathlib import Path
from urllib.request import urlopen,Request
from concurrent.futures import ThreadPoolExecutor
files=[]
for n in ('a370','elgordo','macs0416','a2744','macs0647'):
 for row in csv.DictReader(open(n+'.csv')):
  obs=row['obs_id'];pgm,obsnum=obs.split('_')[0].split('-')
  if obs.endswith(('f356w','f444w','f200w','f150w')): files.append((n,obs,float(row['t_min']),'https://stpubdata.s3.amazonaws.com/jwst/public/'+pgm+'/L3/t/'+obsnum+'/'+obs+'_i2d.fits'))
def probe(item):
 n,obs,m,u=item
 try:
  with urlopen(Request(u,headers={'Range':'bytes=0-2879'}),timeout=12) as r: status=r.status;s=r.headers.get('Content-Range')
 except Exception as e: status=getattr(e,'code',0); s=''
 return [n,obs,m,u,status,s]
rows=list(ThreadPoolExecutor(18).map(probe,files));w=csv.writer(open('cluster_product_access.csv','w'));w.writerow(['field','product','mjd','url','http','range']);w.writerows(rows)
print([(n,obs,round(m)) for n,obs,m,u,st,s in rows if st==206][:70])
