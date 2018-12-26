import pytest

import numpy as np
import astropy.units as u

from astropy.modeling.powerlaws import PowerLaw1D

from pyxsis.binspectrum import XBinSpectrum, group_channels
from pyxsis.bkgspectrum import XBkgSpectrum

def test_init():
    test = XBkgSpectrum([]*u.angstrom, []*u.angstrom,
                        []*u.ct, 1.0*u.second)
    return

@pytest.mark.parametrize('bscal', [0.5, 'array'])
def test_custom_bkg(bscal):
    bin_edges = np.arange(1.0, 10.0, 0.1) * u.angstrom
    if bscal == 'array':
        backscale = np.random.rand(len(bin_edges)-1) + 0.1
    else:
        backscale = bscal

    # Make a fake power law spectrum
    amp0, alpha0 = 3.e-4, 1.0
    powlaw0 = PowerLaw1D(amplitude=amp0, alpha=alpha0, x_0=1.e3)

    bin_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts = np.random.poisson(lam=powlaw0(bin_mid.value), size=len(bin_mid)) * u.ct

    # Initialize the spectrum
    test_bkg = XBkgSpectrum(bin_edges[:-1], bin_edges[1:], counts,
                            exposure=1.0*u.second, backscale=backscale)
    return test_bkg

def test_assign_bkg(bscal='array'):
    test_bkg = test_custom_bkg(bscal)

    # Make a fake power law spectrum for source
    amp0, alpha0 = 3.e-3, -2.0
    powlaw0 = PowerLaw1D(amplitude=amp0, alpha=alpha0, x_0=1.e3)
    bin_mid = test_bkg.spectral_axis.to(u.keV, equivalencies=u.spectral())
    counts  = np.random.poisson(lam=powlaw0(bin_mid.value), size=len(bin_mid)) * u.ct

    test_spec = XBinSpectrum(test_bkg.bin_lo, test_bkg.bin_hi, counts,
                             exposure=1.0*u.second)
    test_spec.assign_bkg(test_bkg)
    return test_spec

# Test binning by channels with and without notice
@pytest.mark.parametrize('use_notice', [True, False])
@pytest.mark.parametrize('bscal', [0.5, 1.0, 'array'])
def test_group_bkg(use_notice, bscal):
    test_spec = test_assign_bkg(bscal)
    if use_notice:
        emin, emax = 3.0 * u.keV, 5.0 * u.keV
        test_spec.notice_range(emin, emax)

    NCHAN = 30
    group_channels(test_spec, NCHAN)
    blo, bhi, cts, cts_err = test_spec.binned_counts()
    bblo, bbhi, bcts, bcts_err = test_spec.binned_bkg()

    # Make sure there are the correct number of counts
    temp_bkg = test_spec.bkg.backscale * test_spec.bkg.counts
    bcounts  = temp_bkg[test_spec.notice]
    np.testing.assert_almost_equal(np.sum(bcounts.value), np.sum(bcts.value), 5)

    assert np.all(blo == bblo)
    assert np.all(bhi == bbhi)

    assert bcts.unit == u.ct
    assert bcts_err.unit == u.ct

    assert len(cts) == len(bcts)
    assert len(cts_err) == len(bcts_err)
