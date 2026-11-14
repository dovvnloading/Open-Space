"""Conditional first-post-Newtonian rate and classical apsidal term, not estimates of mass."""
import json
import numpy as np
from geometry import HERE,P
# Gimenez (1985), eq 2, units degrees per anomalistic cycle, masses solar, P days:
# https://adsabs.harvard.edu/pdf/1985ApJ...297..405G
# Classical equation (rotation aligned) in Claret & Gimenez (1999), eqs. 1–4.
def main():
    out=[]
    allfits=json.loads((HERE/'dilution.json').read_text())['1']
    for row in allfits:
        if row['third']>.301:continue
        _,_,U,e,_,rsum,k,*_=row['pars']
        # For prograde rotation 1/P_anom=1/P_sid - 1/U.
        U_days=U*365.25
        p_anom=1/(1/P-1/U_days)
        w=360*p_anom/U_days
        gr_coeff=5.447e-4/((1-e**2)*p_anom**(2/3))
        required_mass=(w/gr_coeff)**1.5
        r1=rsum/(1+k);r2=k*r1
        mgrid=[]
        for total in (5.,10.,20.,50.,100.):
            gr=gr_coeff*total**(2/3)
            f=(1-e*e)**-2
            g=(1+1.5*e*e+.125*e**4)/(1-e*e)**5
            # Equal masses, spin rates synchronous at mean motion, and equal effective k2.
            csum=(r1**5+r2**5)*(2*f+15*g)
            mgrid.append(dict(total_solar_mass=total,relativistic_fraction=gr/w,
                            implied_effective_k2=(w-gr)/(360*csum)))
        out.append(dict(assumed_dilution=row['third'],p_anomalistic_days=p_anom,
                        apsidal_degrees_per_cycle=w,gr_only_total_solar_mass=required_mass,
                        assumed_mass_examples=mgrid))
    (HERE/'rate.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
