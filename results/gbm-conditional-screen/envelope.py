"""Shared-envelope fit and conditional nuisance bootstrap; engineering diagnostic only."""
import numpy as np
from simulator import DT,N,score,arrival_counts

def fit(counts,width):
 x=np.asarray(counts,float)
 if x.ndim!=2 or x.shape[0]!=2 or x.shape[1]<N or width not in (9,33,65) or np.any(x<0):raise ValueError('invalid fit input')
 weights=np.maximum(x.sum(axis=1),1);scale=weights/weights.sum();shared=x.sum(axis=0);kernel=np.ones(width)/width;pad=width//2;envelope=np.convolve(np.pad(shared,(pad,pad),mode='edge'),kernel,mode='valid');envelope=np.maximum(envelope,1e-3)
 # Fit common AR(1) residual from jointly observed, smoothed rates; this is an
 # approximate nuisance parameter and is not a confidence-calibrated estimator.
 resid=(shared-envelope)/np.sqrt(envelope);phi=float(np.clip(np.dot(resid[1:],resid[:-1])/max(np.dot(resid[:-1],resid[:-1]),1e-12),-.95,.95));sigma=float(np.clip(np.std(resid[1:]-phi*resid[:-1])*.05,.005,.3));return {'expected':np.maximum(scale[:,None]*envelope,1e-3),'phi':phi,'sigma':sigma,'width':width}

def bootstrap(rng,counts,n,width):
 """Refit envelope in every surrogate; bootstrap result conditional on fitted family."""
 x=np.asarray(counts,float);fit0=fit(x,width);L=x.shape[1]
 def max_score(y):
  expected=fit(y,width)['expected'];w=L//N
  if w<1:raise ValueError('insufficient window')
  return score(np.stack([y[:,i*N:(i+1)*N] for i in range(w)]),np.stack([expected[:,i*N:(i+1)*N] for i in range(w)]))
 observed=max_score(x);stats=[]
 for _ in range(n):
  latent=np.zeros(L);phi=float(np.clip(fit0['phi']+rng.uniform(-.1,.1),-.95,.95));sigma=fit0['sigma']*rng.uniform(.5,1.5)
  for j in range(1,L):latent[j]=phi*latent[j-1]+rng.normal(0,sigma)
  shared=np.exp(latent);sim=np.stack([arrival_counts(rng,fit0['expected'][i]*shared/DT) for i in range(2)]);stats.append(max_score(sim))
 return {'observed':observed,'simulated_maxima':stats,'monte_carlo_p':(1+int(np.sum(np.array(stats)>=observed)))/(n+1),'width':width,'phi':fit0['phi'],'sigma':fit0['sigma'],'scope':'conditional toy model only; model choice and nuisance uncertainty not calibrated'}

def joint_bootstrap(rng,bursts,n,width):
 """Resample each burst's common latent process; one joint max across both bursts."""
 fits=[fit(x,width) for x in bursts]
 def measure(x):
  f=fit(x,width);L=x.shape[1];w=L//N
  return score(np.stack([x[:,i*N:(i+1)*N] for i in range(w)]),np.stack([f['expected'][:,i*N:(i+1)*N] for i in range(w)]))
 observed=max(measure(x) for x in bursts);maxima=[]
 for _ in range(n):
  vals=[]
  for x,f in zip(bursts,fits):
   L=x.shape[1];z=np.zeros(L);phi=float(np.clip(f['phi']+rng.uniform(-.1,.1),-.95,.95));sigma=f['sigma']*rng.uniform(.5,1.5)
   for j in range(1,L):z[j]=phi*z[j-1]+rng.normal(0,sigma)
   common=np.exp(z);y=np.stack([arrival_counts(rng,f['expected'][i]*common/DT) for i in range(2)]);vals.append(measure(y))
  maxima.append(max(vals))
 return {'observed':observed,'maxima':maxima,'p':(1+int(np.sum(np.array(maxima)>=observed)))/(n+1),'width':width,'scope':'conditional nuisance family; no calibrated physical significance'}
