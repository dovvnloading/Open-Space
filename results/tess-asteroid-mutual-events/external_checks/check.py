"""Test the two Sector 3 templates on MPC and independent optical measurements."""
import csv
import gzip
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent


def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle))


def write(name, fields, rows):
    with (HERE/name).open('w') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(fields.split())
        writer.writerows(rows)


def field(rows, name):
    return np.array([float(row[name]) for row in rows])


def for_object(rows, number):
    return [row for row in rows if row['mpc'] == str(number)]


def harmonics(time, period):
    angle = time*(2*np.pi*24/period)
    return np.stack([func(k*angle) for k in range(1,6) for func in (np.sin,np.cos)],axis=1)


def later_visits(rows, spin, repeat, phase, value):
    rows = [row for row in rows if row['station']=='C57']
    all_t = field(rows,'jd')
    all_y = field(rows,'mag')
    days = np.unique(all_t.astype(int))
    output = []
    for group in np.split(days,np.where(np.diff(days)>10)[0]+1):
        select = np.isin(all_t.astype(int),group)
        if select.sum()<100:
            continue
        start = min(all_t[select])
        t = all_t[select]-start
        mag = all_y[select]
        keep = np.ones(len(t),dtype=bool)
        def residual(period):
            basis = np.c_[np.ones(len(t)),t/30,(t/30)**2,harmonics(t,period)]
            coeff = np.linalg.lstsq(basis[keep],mag[keep],rcond=None)[0]
            return mag-basis@coeff
        for _ in range(4):
            grid = np.linspace(spin*.995,spin*1.005,161)
            best = min(grid,key=lambda p: np.sum(residual(p)[keep]**2))
            resid = residual(best)
            median = np.median(resid[keep])
            scatter = 1.4826*np.median(abs(resid[keep]-median))
            keep = abs(resid-median)<max(.4,5*scatter)
        t,resid = t[keep],resid[keep]
        trial_phases = np.arange(100)/100
        template = np.stack([np.interp((24*t/repeat-p+.5)%1-.5,phase,value,period=1)
                             for p in trial_phases])
        def amplitudes(select):
            model = template[:,select]
            model = model-model.mean(1)[:,None]
            y = resid[select]
            return model@(y-y.mean())/np.sum(model*model,axis=1)
        half = t<np.median(t)
        first,second = amplitudes(half),amplitudes(~half)
        index = np.argmax(first)
        output.append([round(start,5),int(select.sum()),len(resid),round(best,6),round(scatter,4),
            round(max(amplitudes(np.ones(len(t),bool))),3),trial_phases[index],
            round(first[index],3),trial_phases[np.argmax(second)],round(max(second),3),
            round(second[index],3)])
    return output


def ground_visits(rows, geom, tess, spin, repeat, epoch, phase, value):
    origin = field(tess,'jd').min()
    end = field(tess,'jd').max()
    def event(t):
        return np.interp(((t-epoch)*24/repeat+.5)%1-.5,phase,value,period=1)
    ts = field(tess,'jd')
    angle = field(tess,'phase_deg')
    basis = np.c_[np.ones(len(tess)),angle,angle**2,harmonics(ts-origin,spin),event(ts)]
    mred = field(tess,'mag')-5*np.log10(field(tess,'r_au')*field(tess,'delta_au'))
    spin_amp = np.linalg.lstsq(basis,mred,rcond=None)[0][3:-1]
    output = []
    for pad in (0,7,14,25):
        g = [r for r in rows if r['station'] not in ('C57','I41') and
             origin-pad<float(r['jd'])<end+pad]
        t = field(g,'jd')
        group = np.array([r['station']+'_'+r['band'] for r in g])
        low = np.array([np.interp(t,field(geom,'jd'),field(geom,col))
                        for col in ('r_au','delta_au','phase_deg')])
        alpha = low[2]
        station_filters = np.stack([group==f for f in sorted(set(group))],axis=1)
        X = np.c_[station_filters,alpha,np.minimum(alpha,7),harmonics(t-origin,spin)@spin_amp,event(t)]
        y = field(g,'mag')-5*np.log10(low[0]*low[1])
        coeff = np.linalg.pinv(X)@y
        err = y-X@coeff
        night = t.astype(int)
        nights = np.unique(night)
        scores = np.stack([X[night==d].T@err[night==d] for d in nights])
        bread = np.linalg.pinv(X.T@X)
        correction = len(nights)/(len(nights)-1)*(len(y)-1)/(len(y)-len(coeff))
        se = np.sqrt(np.diag(bread@scores.T@scores@bread)*correction)
        output.append([pad,len(y),len(nights),round(coeff[-2],4),round(se[-2],4),
                       round(coeff[-1],4),round(se[-1],4),round(np.std(err),4)])
    return output


def main():
    obs = read(HERE/'observations.csv.gz')
    geom = read(HERE/'geometry.csv.gz')
    source = read(HERE/'tess_spin_input.csv.gz')
    fold = read(PARENT/'fold.csv')
    measured = read(PARENT/'measurements.csv')
    c57,ground = [],[]
    for number in (29489,3177):
        d = for_object(measured,number)[0]
        phase_data = for_object(fold,number)
        phase,value = field(phase_data,'phase'),field(phase_data,'median_mag')
        spin,repeat,epoch = [float(d[field_name]) for field_name in
                             ('spin_h','recurrence_h','first_center_jd')]
        rows = for_object(obs,number)
        for result in later_visits(rows,spin,repeat,phase,value):
            c57.append([number]+result)
        for result in ground_visits(rows,for_object(geom,number),for_object(source,number),
                                    spin,repeat,epoch,phase,value):
            ground.append([number]+result)
    write('later_tess.csv','mpc jd_start n_input n_used spin_h scatter_mag '
          'full_max_scale first_phase first_scale second_best_phase second_best_scale '
          'second_scale_fixed_first_phase',c57)
    write('ground.csv','mpc pad_days n n_nights spin_scale spin_se event_scale event_se '
          'rms_mag',ground)


if __name__=='__main__':
    main()
