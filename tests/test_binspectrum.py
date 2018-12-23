import pytest

import numpy as np
import astropy.units as u

from astropy.modeling.powerlaws import PowerLaw1D

from pyxsis.binspectrum import XBinSpectrum, group_channels, group_mincounts

def test_init():
    test = XBinSpectrum([]*u.angstrom, []*u.angstrom,
                        []*u.ct, 1.0*u.second)
    return

def test_custom_spec():
    bin_edges = np.arange(1.0, 10.0, 0.1) * u.angstrom

    # Make a fake power law spectrum
    amp0, alpha0 = 3.e-3, -2.0
    powlaw0 = PowerLaw1D(amplitude=amp0, alpha=alpha0, x_0=1.e3)

    bin_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts = np.random.poisson(lam=powlaw0(bin_mid.value), size=len(bin_mid)) * u.ct

    # Initialize the spectrum
    test_spec = XBinSpectrum(bin_edges[:-1], bin_edges[1:], counts,
                             exposure=1.0*u.second)
    return test_spec

def test_notice():
    test_spec = test_custom_spec()

    # test that the noticed values conform to the intervals (keV)
    emin, emax = 3.0 * u.keV, 5.0 * u.keV
    test_spec.notice_range(emin, emax)
    test_notice_keV = test_spec.spectral_axis[test_spec.notice].to(u.keV, equivalencies=u.spectral())
    assert np.all(test_notice_keV >= emin)
    assert np.all(test_notice_keV <= emax)

    # Test that bin edges are within notice range
    bunit = test_spec.bin_lo.unit
    nlo   = test_spec.bin_lo[test_spec.notice]
    nhi   = test_spec.bin_lo[test_spec.notice]
    assert np.all(nlo >= min(emin.to(bunit, equivalencies=u.spectral()),
                             emax.to(bunit, equivalencies=u.spectral())))
    assert np.all(nhi <= max(emin.to(bunit, equivalencies=u.spectral()),
                             emax.to(bunit, equivalencies=u.spectral())))

    # test that the noticed values conform to the intervals (Angstroms)
    amax = emin.to(u.angstrom, equivalencies=u.spectral())
    amin = emax.to(u.angstrom, equivalencies=u.spectral())
    test_spec.notice_range(amin, amax)
    test_notice_angs = test_spec.spectral_axis[test_spec.notice].to(u.angstrom, equivalencies=u.spectral())
    assert np.all(test_notice_angs >= amin)
    assert np.all(test_notice_angs <= amax)

    bunit = test_spec.bin_lo.unit
    nlo   = test_spec.bin_lo[test_spec.notice]
    nhi   = test_spec.bin_lo[test_spec.notice]
    assert np.all(nlo >= min(amin.to(bunit, equivalencies=u.spectral()),
                             amax.to(bunit, equivalencies=u.spectral())))
    assert np.all(nhi <= max(amin.to(bunit, equivalencies=u.spectral()),
                             amax.to(bunit, equivalencies=u.spectral())))


# Test binning by channels with and without notice
@pytest.mark.parametrize('use_notice', [True, False])
@pytest.mark.parametrize('nchan', [10, 30, 10000])
def test_group_channels(use_notice, nchan):
    test_spec = test_custom_spec()
    if use_notice:
        emin, emax = 3.0 * u.keV, 5.0 * u.keV
        test_spec.notice_range(emin, emax)

    group_channels(test_spec, nchan)
    blo, bhi, cts, cts_err = test_spec.binned_counts()

    notice  = test_spec.notice
    binning = test_spec.binning[notice]
    counts  = test_spec.counts[notice].value
    bin_lo  = test_spec.bin_lo[notice].value
    bin_hi  = test_spec.bin_hi[notice].value

    assert len(blo) == (max(binning) - min(binning) + 1)
    assert len(bhi) == (max(binning) - min(binning) + 1)
    assert len(cts) == (max(binning) - min(binning) + 1)
    assert np.sum(cts.value) == np.sum(counts)  # Make sure no counts are lost
    assert np.all(blo - bhi < 0.0)

# Test binning by counts with and without notice
@pytest.mark.parametrize('use_notice', [True, False])
@pytest.mark.parametrize('mcounts', [1, 10, 30])
def test_group_channels(use_notice, mcounts):
    test_spec = test_custom_spec()
    if use_notice:
        emin, emax = 3.0 * u.keV, 5.0 * u.keV
        test_spec.notice_range(emin, emax)

    group_mincounts(test_spec, mcounts)
    blo, bhi, cts, cts_err = test_spec.binned_counts()

    notice  = test_spec.notice
    binning = test_spec.binning[notice]
    counts  = test_spec.counts[notice].value
    bin_lo  = test_spec.bin_lo[notice].value
    bin_hi  = test_spec.bin_hi[notice].value

    assert len(blo) == (max(binning) - min(binning) + 1)
    assert len(bhi) == (max(binning) - min(binning) + 1)
    assert len(cts) == (max(binning) - min(binning) + 1)
    assert np.sum(cts.value) == np.sum(counts)  # Make sure no counts are lost
    assert np.all(blo - bhi < 0.0)
