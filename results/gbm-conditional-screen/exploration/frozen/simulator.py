"""Conditional joint-detector surrogate simulation; no real data access."""
import numpy as np
DT=1/64;N=256;FREQ=np.fft.rfftfreq(N,DT);F=(FREQ>=5)&(FREQ<=32)
def score(counts,mean):
 counts=np.asarray(counts,float);mean=np.asarray(mean,float)
 if counts.shape!=mean.shape or counts.ndim!=3 or counts.shape[1:]!=(2,N) or (mean<=0).any():raise ValueError('shape/mean')
 x=counts-mean;p=2*np.abs(np.fft.rfft(x,axis=2))**2/mean.sum(axis=2,keepdims=True);p[:,:,-1]*=.5
 joint=np.minimum(p[:,0,F],p[:,1,F]);return float(np.max(joint))
def arrival_counts(rng,rates,deadtime=2.6e-6,overflow_prob=0.):
 """Within-bin Poisson events and nonparalyzable per-detector loss; bin boundaries share deadtime."""
 out=np.zeros(len(rates),int);last=-1e100
 for i,rate in enumerate(rates):
  n=rng.poisson(rate*DT);times=i*DT+rng.random(n)*DT;times.sort()
  for t in times:
   if t>=last:
    out[i]+=1;last=t+(1e-5 if rng.random()<overflow_prob else deadtime)
 return out
def bin_modulation(t,frequency,phase=0.):
 return np.sinc(frequency*DT)*np.sin(2*np.pi*frequency*(np.asarray(t)+DT/2)+phase)

def simulate(rng,kind='smooth',injection=0.,frequency=10.,phase=0.,nwindows=7,deadtime=2.6e-6):
 """Shared red-noise realization in each burst window; detector counting noise separate."""
 x=np.arange(N)*DT;means=[];draw=[]
 for w in range(nwindows):
  phase=(x-2)/(.9 if kind=='smooth' else .25);env=500+2500*np.exp(-.5*phase**2)
  if kind=='pulse_train':env+=sum(300*np.exp(-.5*((x-q)/.07)**2) for q in (.7,1.3,2.,2.7,3.3))
  latent=np.zeros(N);latent[0]=previous if w not in (0,2) else rng.normal(0,.08)
  for k in range(1,N):latent[k]=.95*latent[k-1]+rng.normal(0,.025)
  previous=latent[-1]
  intrinsic=np.clip(env*np.exp(latent)*(1+injection*bin_modulation(x,frequency,phase)),1,None)
  m=np.stack([env,env*.7])*DT;means.append(m);draw.append(np.stack([arrival_counts(rng,intrinsic,deadtime),arrival_counts(rng,intrinsic*.7,deadtime)]))
 return np.asarray(draw),np.asarray(means)
def monte_carlo(seed,n,kind='smooth',injection=0.):
 rng=np.random.default_rng(seed);return np.array([score(*simulate(rng,kind=kind,injection=injection)) for _ in range(n)])
