import numpy as np
import astropy.units as u
from astropy.modeling.powerlaws import PowerLaw1D

from pyxsis.binspectrum import XBinSpectrum

def test_init():
    test = XBinSpectrum([]*u.angstrom, []*u.angstrom,
                        []*u.ct, 1.0*u.second)
    return

def test_generic():
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
