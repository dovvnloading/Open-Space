"""Illustrative bolometric fraction a blackbody emits between 0.75 and 5 microns."""
import csv
import math
from pathlib import Path
from scipy.integrate import quad

# Central estimates in Hephaistos III tables 1 and 3, https://arxiv.org/abs/2607.25701.
# These reproduce published hypotheses including cases now identified as background galaxies.
TEMPERATURES=[('A',138,6),('B',275,40),('C',187,16),('D',178,20),('E',180,26),
              ('F',176,21),('G',100,9),('H',130,21),('I',99,20),('J',114,28)]

def fraction(k):
    a=14387.76877/(5*k); b=14387.76877/(.75*k)
    return 15/math.pi**4 * quad(lambda x:x**3*math.exp(-x)/(-math.expm1(-x)),a,b,epsabs=1e-20)[0]

if __name__=='__main__':
    with (Path(__file__).parent/'thermal_leverage.csv').open('w') as f:
        w=csv.writer(f,lineterminator='\n')
        w.writerow(['label','published_T_K','published_T_err_K','bolometric_fraction_0p75_to_5um','fraction_at_T_minus_err','fraction_at_T_plus_err'])
        for name,t,err in TEMPERATURES:
            w.writerow([name,t,err,*[f'{fraction(value):.8g}' for value in [t,t-err,t+err]]])
