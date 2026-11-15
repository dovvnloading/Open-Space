"""Independent one-cycle holdout of the 2018 edge using the next XMM visit."""
from analyze import *
EDGE=58252.13450
P=(LATE_EGRESS-EPOCH)/507
def one(f):
    prefix,p,g=read_visit(f)
    predicted=EDGE+2*P
    out=[]
    # Adjacent pre/post bins exclude the complete first-edge profile interval.
    for label,l,r in [('before',predicted-.28,predicted-.025),('after',predicted+.025,predicted+.24)]:
        inside=np.maximum(0,np.minimum(r,g[:,1])-np.maximum(l,g[:,0]))
        dur=inside.sum()*86.4
        if dur<1:continue
        m=p.energy_ev.between(900,7000)&p.mjd_tdb_geo.between(l,r)&(abs(p.x_arcsec-NEIGHBOR[0]/2)<19)&(abs(p.y_arcsec-NEIGHBOR[1]/2)<19)
        v=fit(p.x_arcsec[m].to_numpy(),p.y_arcsec[m].to_numpy())
        if v is not None:
            out.append(dict(detector=prefix[11:13],predicted_mjd=predicted,window=label,window_start=l,window_end=r,exposure_ks=dur,ulx1_rate=v[0]/dur,ulx1_error=v[1]/dur,ulx2_rate=v[2]/dur,delta_c_source=v[3]))
    return out
if __name__=='__main__':
    result=[r for f in sorted((HERE/'inputs').glob('P0830191401*-photons.csv.gz')) for r in one(f)]
    print(pd.DataFrame(result).round(3).to_string(index=False))
    pd.DataFrame(result).to_csv(HERE/'xmm-2018-recurrence.csv',index=False,float_format='%.8f')
