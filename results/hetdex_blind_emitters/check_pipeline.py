"""Small deterministic regression check of the actual selection entry points."""
import numpy as np

from scan import screen_chunk
from spectral import simulate, search_batch, shuffle_fibres


def main():
    wave = np.arange(3470., 5541., 2.)
    rng = np.random.default_rng(5727)
    error = np.ones((4, len(wave)))
    flux = rng.normal(size=error.shape)
    centers = (4511., 4553., 4963., 5460.)
    for row, center in zip(flux, centers):
        row += 12 * np.exp(-0.5*((wave-center)/2.6)**2)
    peaks = search_batch(wave, flux, error, density_unit=1e-17)
    for idx, center in enumerate(centers[:3]):
        assert any(i==idx and abs(p["wave"]-center)<3 for i,p in peaks), (idx, peaks)
    assert not any(i==3 and abs(p["wave"]-5460)<3 for i,p in peaks)
    mock, mock_err=shuffle_fibres(flux,error,seed=12)
    assert np.array_equal(np.sort(mock,axis=0),np.sort(flux,axis=0))
    for y,e in ((flux,error),(mock,mock_err)):
        counts,_=screen_chunk(wave,y,e,np.full(4,2.2),0,'fixture')
        assert counts['spectra']==4
    assert len(simulate(wave,flux[0],error[0],2,n_inject=6,density_unit=1e-17)['injections'])==6
    print('Spectral recovery, calibration masking, shuffled selection and injection entry points passed.')


if __name__=='__main__':
    main()
