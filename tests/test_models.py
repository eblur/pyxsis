import numpy as np
import astropy.units as u
from xpysis.models import *

def test_model_superclass():
    names = ['a','b','c']
    vals  = [1.0, 2.0, 3.0]
    lims  = [(0,10), (0,10), (0,10)]
    units = ['s', 'cm^2', 'keV']
    m = Model(names, vals, lims, units)
    for k in names:
        assert isinstance(m[k], u.Quantity)
