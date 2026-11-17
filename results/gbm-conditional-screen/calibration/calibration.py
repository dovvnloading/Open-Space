"""Synthetic-only retrospective pulse-family calibration. No input data access."""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from simulator import DT,N,arrival_counts,score,bin_modulation
from envelope import fit

def draw(seed,variant='training',inject=0.):
 rng=np.random.default_rng(seed);pieces=[];models=[]
 for windows in (2,5):
  length=windows*N;t=np.arange(length)*DT;rate=np.ones(length)*rng.uniform(350,700)
  for _ in range(rng.integers(8,20)):
   center=rng.uniform(t[0],t[-1]);width=rng.uniform(.04,.18) if variant!='narrow' else rng.uniform(.012,.05);amp=rng.uniform(400,2200);rate+=amp*np.exp(-.5*((t-center)/width)**2)
  red=np.zeros(length)
  for j in range(1,length):red[j]=.95*red[j-1]+rng.normal(0,.025)
  common=np.maximum(rate*np.exp(red)*(1+inject*bin_modulation(t,10.,.3)),1)
  # Detector-dependent response includes a fixed relative scale and small slow drift.
  r1=common;r2=common*(.65+.05*np.sin(2*np.pi*t/13))
  x=np.stack([arrival_counts(rng,r1),arrival_counts(rng,r2)]);e=fit(x,9)['expected']
  pieces.extend([x[:,i*N:(i+1)*N] for i in range(windows)]);models.extend([e[:,i*N:(i+1)*N] for i in range(windows)])
 return score(np.stack(pieces),np.stack(models))

def summarize(n=100):
 train=np.array([draw(10000+i) for i in range(n)]);test=np.array([draw(20000+i) for i in range(n)]);alt=np.array([draw(30000+i,'narrow') for i in range(n)]);sig=np.array([draw(40000+i,inject=.2) for i in range(n)]);threshold=float(np.quantile(train,.95,method='higher'))
 return {'n_each':n,'threshold_training_95pct':threshold,'heldout_same_family_crossings':int((test>=threshold).sum()),'heldout_narrow_family_crossings':int((alt>=threshold).sum()),'heldout_injected_crossings':int((sig>=threshold).sum()),'scope':'synthetic family only; not a calibrated real-burst p-value'}
