"""Explicitly authorized joint detector phase; no automatic execution or network."""
import argparse,json
from pathlib import Path
import numpy as np
from astropy.io import fits
from manifest import RAW,sha,verify
from fits_reader import metadata,decode

SPECS=[('261002867','n5','NAI_05',812666888.812348,-.12803,2),('261002867','n2','NAI_02',812666888.812348,-.12803,2),('261007236','nb','NAI_11',813044394.710608,2.84802,5),('261007236','n8','NAI_08',813044394.710608,2.84802,5)]
def authorize(path,digest,manifest_digest):
 if sha(path)!=digest:raise ValueError('authorization digest')
 a=json.loads(Path(path).read_text())
 if a.get('authorized') is not True or a.get('phase')!='joint_tte' or a.get('manifest_sha256')!=manifest_digest:raise ValueError('authorization content')
def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--manifest-sha256',required=True);p.add_argument('--authorization',required=True);p.add_argument('--authorization-sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();m=verify(a.manifest,a.manifest_sha256);authorize(a.authorization,a.authorization_sha256,a.manifest_sha256)
 if a.output.exists() or a.output.resolve().drive.upper()!='D:':raise ValueError('output collision or not D drive')
 handles=[]
 try:
  for bid,det,name,trig,start,n in SPECS:
   filename=f'glg_tte_{det}_bn{bid}_v00.fit'
   if filename not in m['tte_sha256']:raise ValueError('input not in manifest')
   h=fits.open(RAW/filename,memmap=True);handles.append(h);metadata(h,name,bid,trig,start,start+n*4)
  # All identities, hashes, GTIs, calibration schemas and authorization are
  # checked before any events table is accessed.
  output=[]
  for h,(bid,det,name,trig,start,n) in zip(handles,SPECS):
   ch,selected=metadata(h,name,bid,trig,start,start+n*4);t,_,_=decode(h,ch,selected,trig,start,start+n*4);counts=np.histogram(t,bins=start+np.arange(n*256+1)/64)[0];output.append({'burst':bid,'detector':det,'counts':counts.tolist(),'start_s':start,'nwindows':n})
  with a.output.open('x',encoding='utf-8') as f:json.dump({'manifest_sha256':a.manifest_sha256,'statistic':'joint maximum; no detector independence','series':output},f,allow_nan=False);f.write('\n')
 finally:
  for h in handles:h.close()
if __name__=='__main__':main()
