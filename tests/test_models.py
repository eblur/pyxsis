import numpy as np
import astropy.units as u
from xpysis.models import Model, PowerLaw

def test_model_superclass():
    names = ['a','b','c']
    vals  = [1.0, 2.0, 3.0]
    lims  = [(0.,10.), (0.,10.), (0.,10.)]
    units = ['s', 'cm^2', 'keV']
    m = Model(names, vals, lims, units)
    # Check that __getitem__ works
    for k in names:
        assert isinstance(m[k], u.Quantity)
    # Check that initializing with an out-of-bound parameter does not work
    # oops, need to look up catching assertion error with py.test
    #assert Model(names, [-10., 2.0, 3.0], lims, units) is None
    # Check that update works in general
    # and that updating to an out-of-bound parameter will set it to limit value
    new_vals = dict(zip(['a','b','c'], [2.5, -10.0, 25.0]))
    m.update(new_vals)
    assert m['a'].value == 2.5
    assert m['b'].value == 0.0
    assert m['c'].value == 10.0
    m.update({'a':3.0})
    assert m.vals['a'] == 3.0
    # Check that placeholder method is there
    assert m.calculate(1.0, 2.0) is None

def test_powerlaw():
    # Check that it initializes
    ener = np.linspace(1.0, 10.0)
    e_lo, e_hi = ener[:-1], ener[1:]
    m = PowerLaw()
    # Check that the flux can be evaluated
    flux = m.calculate(e_lo, e_hi)
    assert flux is not None
    assert len(flux) == len(e_lo)
    assert flux.unit == u.Unit(m.units['norm'])
    # Check that doubling the norm parameter doubles the Flux
    new_norm = m['norm'].value * 2.0
    m2 = PowerLaw(new_norm)
    f2 = m2.calculate(e_lo, e_hi)
    assert np.sum(f2) == 2.0 * np.sum(flux)
    # Check that updating parameters works
    m.update({'norm':new_norm})
    assert np.sum(m.calculate(e_lo, e_hi)) == np.sum(f2)
