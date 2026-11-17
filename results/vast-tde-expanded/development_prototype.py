"""Review prototype: development-only calibration and frozen prediction.
No archive reader is provided; caller must supply development-only rows.
"""
import numpy as np

def prepare(source_id,packet,rows):
 cfg=packet['sources'][str(source_id)];allowed={i for group in cfg['development'] for i in group};held={i for group in cfg['holdout'] for i in group}
 if any(int(r['image_id']) in held or int(r['image_id']) not in allowed for r in rows):raise ValueError('non-development image encountered')
 if len({(r['source'],r['image_id']) for r in rows})!=len(rows):raise ValueError('duplicate source-image')
 return rows

def constant_prediction(values,errors):
 x=np.asarray(values,float);s=np.asarray(errors,float)
 if len(x)<4 or x.shape!=s.shape or not np.isfinite(x).all() or not np.isfinite(s).all() or np.any(s<=0):raise ValueError('insufficient valid development sessions')
 w=1/s**2;mean=float(np.sum(w*x)/sum(w));chi=float(np.sum(w*(x-mean)**2));inflate=max(1.,chi/(len(x)-1));return dict(mean=mean,variance=float(inflate/sum(w)),chi2=chi,dof=len(x)-1)

def conditional_gate(pred,holdout,holdout_cov):
 h=np.asarray(holdout,float);c=np.asarray(holdout_cov,float)
 if h.shape!=(2,) or c.shape!=(2,2) or not np.isfinite(h).all() or not np.isfinite(c).all() or not np.allclose(c,c.T) or np.linalg.eigvalsh(c).min()<=0:raise ValueError('invalid heldout covariance')
 total=c+pred['variance']*np.ones((2,2));z=(h-pred['mean'])/np.sqrt(np.diag(total));return {'z':z.tolist(),'conditional_rise':bool(np.all(z>=3))}

if __name__=='__main__':
 packet={'sources':{'1':{'development':[[1],[2],[3],[4]],'holdout':[[5],[6]]}}}
 assert len(prepare(1,packet,[{'source':1,'image_id':1}]))==1
 try:prepare(1,packet,[{'source':1,'image_id':5}])
 except ValueError:pass
 else:raise AssertionError('holdout leaked')
 pred=constant_prediction([1,1,1,1],[.1]*4)
 assert not conditional_gate(pred,[1,1],np.diag([.01,.01]))['conditional_rise']
 assert conditional_gate(pred,[1.5,1.5],np.diag([.01,.01]))['conditional_rise']
 print('development guard and synthetic gate passed; no real measurements read')
