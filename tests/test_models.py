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

def test_powerlaw():
    ener = np.linspace(1.0, 10.0)
    e_lo, e_hi = ener[:-1], ener[1:]
    m = PowerLaw()
    flux = m.calculate(e_lo, e_hi)
    assert flux is not None
